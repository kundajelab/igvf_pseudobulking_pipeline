import logging
from concurrent.futures import (
    ThreadPoolExecutor,
    as_completed,
)
from multiprocessing import cpu_count
from pathlib import Path
from threading import Lock

import pandas as pd

from pseudobulk import utils


def _load_and_combine_accession_qc(
    accession: str,
    atac_qc: pd.DataFrame,
    rna_qc_dir: Path,
    output_dir: Path,
    log_lock: Lock,
) -> None:
    """Load and combine ATAC and RNA QC for a given accession.

    Args:
        accession: Analysis accession ID to process.
        atac_qc: DataFrame with all ATAC QC data
        rna_qc_dir: Path to folder with RNA QC TSV.
        logger: Logger instance for logging.
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(name=f"{__package__} combine-accession-qc")
    with log_lock:
        logger.info(f"Processing analysis accession: {accession}")
    # Load combined atac qc (generated per analysis accession, not by pseudobulk)
    # NOTE: don't need raw columns after pseudobulk processing
    atac_qc = atac_qc.loc[atac_qc["analysis_set_accession"] == accession, :].copy()
    atac_qc["found_in_atac"] = True
    with log_lock:
        logger.info(f"Got ATAC QC for {accession} with shape {atac_qc.shape}")
    rna_qc_tsv = rna_qc_dir / f"{accession}.scRNA_all_cells_QC_metrics.tsv"
    if rna_qc_tsv.exists():
        rna_qc = utils.read_csv(rna_qc_tsv)
        rna_qc["found_in_rna"] = True
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
    logger = logging.getLogger(name=f"{__package__} combine-accession-qc")
    num_workers = num_workers if num_workers > 0 else cpu_count()
    logger.info(f"Executing with {num_workers} worker threads.")
    # Load metadata
    logger.info("Loading metadata")
    metadata_df = utils.load_metadata(metadata_loc)
    log_lock = Lock()

    atac_qc = utils.load_atac_qc(
        atac_qc_dir=atac_qc_dir,
        need_tss=False,
        logger=logger,
        log_lock=log_lock,
        usecols=lambda _col: not f"{_col}".startswith("raw-"),
    )

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Combine QC per analysis accession
        futures = [
            executor.submit(
                _load_and_combine_accession_qc,
                accession=accession,
                atac_qc=atac_qc,
                rna_qc_dir=rna_qc_dir,
                output_dir=output_dir,
                log_lock=log_lock,
            )
            for accession in metadata_df["analysis_set_accession"].unique().tolist()
        ]
        for future in as_completed(futures):
            future.result()  # Raise any exceptions that occurred in the worker processes
