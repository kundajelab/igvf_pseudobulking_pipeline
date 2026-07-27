import logging
from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed,
)
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import Lock, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd

from pseudobulk import utils

_LOGGER_NAME = f"{__package__} combine-accession-qc"

# The log lock is held once per worker process rather than passed with every task.
# NOTE: pandas is not safe to use from several threads with the GIL disabled, so the workers are
# processes, which also means the lock has to be a multiprocessing one.
_WORKER_LOG_LOCK: dict[str, utils.LogLock] = {}


def _init_worker(log_lock: utils.LogLock, log_level: int) -> None:
    """Initialize a worker process with logging and the shared log lock."""
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    _WORKER_LOG_LOCK["log_lock"] = log_lock


def _load_and_combine_accession_qc(
    accession: str,
    atac_qc_dir: Path,
    rna_qc_dir: Path,
    output_dir: Path,
) -> None:
    """Load and combine ATAC and RNA QC for a given accession.

    Args:
        accession: Analysis accession ID to process.
        atac_qc_dir: Path to folder with ATAC QC TSVs.
        rna_qc_dir: Path to folder with RNA QC TSV.
        output_dir: Path to folder to write the combined per-accession QC TSV to.
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    log_lock = _WORKER_LOG_LOCK["log_lock"]
    with log_lock:
        logger.info(f"Processing analysis accession: {accession}")
    # Load the atac qc of this accession only (it is generated per analysis accession, not by
    # pseudobulk), so that each worker holds just the QC it needs.
    # NOTE: don't need raw columns after pseudobulk processing
    atac_qc = utils.load_atac_qc(
        atac_qc_dir=atac_qc_dir,
        need_tss=False,
        logger=logger,
        log_lock=log_lock,
        usecols=lambda _col: not f"{_col}".startswith("raw-"),
        accessions={accession},
    )
    atac_qc["found_in_atac"] = True
    with log_lock:
        logger.info(f"Got ATAC QC for {accession} with shape {atac_qc.shape}")
    rna_qc_tsv = rna_qc_dir / f"{accession}.scRNA_all_cells_QC_metrics.tsv"
    if rna_qc_tsv.exists():
        rna_qc = utils.read_csv(rna_qc_tsv)
        rna_qc["found_in_rna"] = True
        rna_qc["rna_read_count"] = rna_qc["rna_read_count"].astype(np.uint64)
        rna_qc["gene_count"] = rna_qc["gene_count"].astype(np.uint64)
        with log_lock:
            logger.info(f"Got RNA QC for {accession} with shape {rna_qc.shape}")
    else:
        rna_qc = pd.DataFrame([], columns=pd.Index(utils.RNA_QC_COLUMNS + ("found_in_rna",)))
    # Combine QC and write out
    combined_qc = utils.merge_rna_and_atac_qc(
        identifier=accession, atac_qc=atac_qc, rna_qc=rna_qc, logger=logger, log_lock=log_lock
    )
    with log_lock:
        logger.info(f"Writing combined QC for {accession}")
    combined_qc.to_csv(f"{output_dir}/{accession}_per_cell_qc.tsv.gz", sep="\t", index=False)


def combine_accession_qc(
    *,
    metadata_loc: Path,
    atac_qc_dir: Path,
    rna_qc_dir: Path,
    output_dir: Path,
    num_workers: int = -1,
):
    """Combine ATAC and RNA QC for each accession.

    Args:
        metadata_loc: Input annotations metadata file path
        atac_qc_dir: Path to folder with ATAC QC TSVs.
        rna_qc_dir: Path to folder with RNA QC TSV.
        output_dir: Path to folder to write per-accession combined QC TSVs.
        num_workers: Number of parallel workers to use. If <=0, use all available cores.
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    num_workers = num_workers if num_workers > 0 else cpu_count()
    logger.info(f"Executing with {num_workers} worker processes.")
    # Load metadata
    logger.info("Loading metadata")
    metadata_df = utils.load_metadata(metadata_loc)

    # note the OOM kill count now, so that a worker that is killed later can be identified as having
    # run out of memory
    oom_kills_before = utils.oom_kill_count()
    try:
        with ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=_init_worker,
            initargs=(Lock(), logger.getEffectiveLevel()),
        ) as executor:
            # Combine QC per analysis accession
            futures = [
                executor.submit(
                    _load_and_combine_accession_qc,
                    accession=accession,
                    atac_qc_dir=atac_qc_dir,
                    rna_qc_dir=rna_qc_dir,
                    output_dir=output_dir,
                )
                for accession in metadata_df["analysis_set_accession"].unique().tolist()
            ]
            for future in as_completed(futures):
                future.result()  # Raise any exceptions that occurred in the worker processes
    except BrokenProcessPool:
        utils.exit_if_oom_killed(executor, oom_kills_before, logger)
        raise
