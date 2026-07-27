import gc
import logging
import tempfile
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed,
)
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import Lock, cpu_count
from pathlib import Path
from typing import Final

import anndata as ad
import numpy as np
import pandas as pd

from pseudobulk import sc_qc, utils
from pseudobulk.types import PseudobulkName
from pseudobulk.utils import load_metadata

_LOGGER_NAME: Final[str] = "pseudobulk-rna"

# Reference data that every task needs, held once per worker process. Passing these as task
# arguments instead would re-pickle them (the metadata is ~300 MB in memory) for every h5ad.
# NOTE: pandas is not safe to use from several threads with the GIL disabled, so the workers are
# processes, and AnnData objects are handed between them as pickle files rather than as arguments.
_WORKER_DATA: dict[str, pd.DataFrame] = {}

# The log lock is held once per worker process rather than passed with every task. It has to be a
# multiprocessing lock, because a threading lock only serializes the threads within one process.
_WORKER_LOG_LOCK: dict[str, utils.LogLock] = {}


def _init_worker(
    gene_ref: pd.DataFrame, metadata_loc: Path, log_lock: utils.LogLock, log_level: int
) -> None:
    """Initialize a worker process with logging and the shared reference data."""
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    _WORKER_DATA["gene_ref"] = gene_ref
    _WORKER_DATA["metadata_df"] = utils.load_metadata(
        metadata_loc=metadata_loc,
        wanted_cols=["analysis_set_accession", "barcode_sample", "pseudobulk_id"],
    )
    _WORKER_LOG_LOCK["log_lock"] = log_lock


def _worker_data(name: str) -> pd.DataFrame:
    """Get reference data held by this worker process."""
    try:
        return _WORKER_DATA[name]
    except KeyError:
        raise RuntimeError(f"worker process has no {name}: it was not initialized") from None


def _worker_log_lock() -> utils.LogLock:
    """Get the log lock held by this worker process."""
    try:
        return _WORKER_LOG_LOCK["log_lock"]
    except KeyError:
        raise RuntimeError("worker process has no log lock: it was not initialized") from None


def _get_h5ad_paths(input: Path, logger: logging.Logger) -> Iterable[Path]:
    """Get paths to raw RNA h5ad files from input.

    Args:
      input: Either a folder, in which case raw input RNA files will be globbed from this folder,
            or a file, with one input RNA file per line.
    Returns:
        Iterable of raw RNA h5ad files.
    """
    if input.is_dir():
        logger.info(f"Finding h5ad files in {input}")
        yield from input.glob("*.h5ad")
    elif input.is_file():
        logger.info(f"Reading input paths from file-of-files: {input}")
        with open(input, "rt") as f_in:
            for h5ad_file in f_in:
                yield Path(h5ad_file)
    else:
        raise ValueError("input is neither a folder nor a file")


def _load_ann_data(h5ad_path: Path, gene_ref: pd.DataFrame) -> ad.AnnData:
    """Load and preprocess h5ad file with raw RNA data."""
    adata = ad.read_h5ad(f"{h5ad_path}", backed="r")
    obs = pd.DataFrame(adata.obs)
    obs["analysis_set_accession"] = h5ad_path.name.split(".", 1)[0]
    obs["analysis_set_accession"] = obs["analysis_set_accession"].astype("category")
    obs["barcode_sample"] = obs.index
    adata.obs = obs
    # Compute QC for each file (analysis_set_accession)
    var = pd.DataFrame(adata.var)
    var.loc[:, "gene_symbol"] = var.index.map(gene_ref["gene_name"]).array
    var.loc[:, "mt"] = var.index.map(gene_ref["mt"]).fillna(False).astype(bool)
    var.loc[:, "ribo"] = var.index.map(gene_ref["ribo"]).fillna(False).astype(bool)
    adata.var = var
    adata.uns["h5ad_path"] = f"{h5ad_path}"
    # TODO: add tests showing equivalence
    # sc.pp.calculate_qc_metrics(
    #     adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    # )
    sc_qc.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"])
    return adata


def _update_and_save_obs(
    adata: ad.AnnData, metadata_df: pd.DataFrame, qc_report_path: Path | None = None
) -> pd.DataFrame:
    """Save QC info for AnnData corresponding to a raw RNA h5ad."""
    obs = pd.DataFrame(adata.obs).rename(
        columns={
            "total_counts": "rna_read_count",
            "n_genes_by_counts": "gene_count",
            "pct_counts_mt": "pct_mito",
            "pct_counts_ribo": "pct_ribo",
        },
    )  # Rename
    obs["rna_read_count"] = obs["rna_read_count"].astype(np.uint64)
    obs["gene_count"] = obs["gene_count"].astype(np.uint64)
    obs["annotated"] = obs["barcode_sample"].isin(set(metadata_df["barcode_sample"]))
    barcodes_to_pseudobulks = utils.map_barcodes_to_pseudobulks(metadata_df)
    obs["pseudobulk_id"] = obs["barcode_sample"].map(
        lambda barcode: barcodes_to_pseudobulks.get(barcode, "null")
    )
    obs["pseudobulk_id"] = obs["pseudobulk_id"].astype("category")
    # restrict to desired columns
    wanted_cols = list(utils.RNA_QC_COLUMNS)
    obs = obs.loc[:, wanted_cols]
    if qc_report_path is not None:
        # Save QC for all cells in analysis accession
        qc_report_path.parent.mkdir(parents=True, exist_ok=True)
        sep = utils.get_sep_from_path(qc_report_path)
        obs.to_csv(f"{qc_report_path}", sep=sep, index=False)
    return obs


def _save_cell_name_to_annotation_mapping(metadata_loc: Path, output_dir: Path) -> None:
    """Get mapping from cell name to annotation and save to TSV."""
    cell_name_to_annotation_df = load_metadata(
        metadata_loc,
        wanted_cols=[
            "pseudobulk_id",
            "cell_name",
            "annotation",
            "CL_id",
            "cell_description",
            "CL_term_name",
            "subsample",
            "cell_qualifier",
        ],
    ).drop_duplicates()
    # note that cell_qualifier may or may not be present, it is an optional field
    cell_name_to_annotation_df.to_csv(
        f"{output_dir}/cell_name_to_annotation_mapping.tsv", sep="\t", index=False
    )


def _save_adata_pseudobulk(
    adata: ad.AnnData,
    pseudobulk_id: PseudobulkName,
    pseudobulk_metadata: pd.DataFrame,
    h5ad_path: Path,
    accession: str,
    temp_pseudobulk_dir: Path,
    log_lock: utils.LogLock,
    logger: logging.Logger,
) -> tuple[PseudobulkName, Path]:
    """Separate out the AnnData for the pseudobulk, and save to a temporary path."""
    with log_lock:
        logger.info(f"\tIn {h5ad_path}: separating out AnnData for pseudobulk_id {pseudobulk_id}")
    barcode_samples = set(pseudobulk_metadata["barcode_sample"])
    wanted_rows = adata.obs_names.isin(barcode_samples)
    adata_path = temp_pseudobulk_dir / f"{pseudobulk_id}.{accession}.h5ad"
    pseudobulk_adata: ad.AnnData = adata[wanted_rows, :].copy(filename=f"{adata_path}")
    # garbage collection is needed because of how AnnData objects are stored.
    del pseudobulk_adata
    gc.collect()
    return PseudobulkName(f"{pseudobulk_id}"), adata_path


def _load_and_qc_h5ad(
    h5ad_path: Path,
    rna_qc_reports_dir: Path,
    temp_pseudobulk_dir: Path,
) -> list[tuple[PseudobulkName, Path]]:
    """Load one raw RNA h5ad, save its QC, and write out one pickle per pseudobulk ID.

    Returns:
        list of (pseudobulk ID, path to the pickled AnnData for that pseudobulk).
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    log_lock = _worker_log_lock()
    metadata_df = _worker_data("metadata_df")
    accession = h5ad_path.name.split(".", 1)[0]
    with log_lock:
        logger.info(f"Processing h5ad file: {h5ad_path} for accession {accession}")
    metadata_df: pd.DataFrame = metadata_df.loc[
        metadata_df["analysis_set_accession"] == accession, :
    ]
    adata = _load_ann_data(h5ad_path, gene_ref=_worker_data("gene_ref"))

    qc_report_path = rna_qc_reports_dir / f"{accession}.scRNA_all_cells_QC_metrics.tsv"
    adata.obs = _update_and_save_obs(adata, metadata_df=metadata_df, qc_report_path=qc_report_path)
    adata.strings_to_categoricals()

    # Get annotated data for each pseudobulk ID
    pseudobulk_adata_paths = [
        _save_adata_pseudobulk(
            adata=adata,
            pseudobulk_id=PseudobulkName(pseudobulk_id),  # ty:ignore[invalid-argument-type]
            pseudobulk_metadata=pseudobulk_metadata,
            h5ad_path=h5ad_path,
            accession=accession,
            temp_pseudobulk_dir=temp_pseudobulk_dir,
            log_lock=log_lock,
            logger=logger,
        )
        # observed=True because pseudobulk_id is a categorical whose categories span every
        # pseudobulk in the run: the default would yield a group for each of them for every
        # accession, almost all empty, and write an h5f5 for each one.
        for pseudobulk_id, pseudobulk_metadata in metadata_df.groupby(
            "pseudobulk_id", sort=False, group_keys=False, as_index=False, observed=True
        )
    ]
    # AnnData objects take part in reference cycles, so nothing but the cycle collector can free
    # them, and it is triggered by the number of allocations rather than by their size. A handful
    # of hundred-MB objects never reaches that threshold, so without collecting here each task
    # leaves its whole h5ad behind, growing with each file without bound, until the task was killed
    # for running out of memory.
    del adata
    gc.collect()
    return pseudobulk_adata_paths


def _load_and_qc_h5ads(
    executor: ProcessPoolExecutor,
    metadata_loc: Path,
    h5ad_paths: Iterable[Path],
    output_dir: Path,
    rna_qc_reports_dir: Path,
    temp_pseudobulk_dir: Path,
) -> defaultdict[PseudobulkName, list[Path]]:
    """Load raw RNA h5ads and separate into pseudobulked AnnData objects.

    Also save observation and QC data for each accession ID into rna_qc_reports_dir.

    Args:
        executor: process pool to load the h5ads with. Its workers must have been initialized with
            _init_worker, which holds the metadata and gene reference data they need.
        metadata_df: metadata DataFrame, as loaded by utils.load_metadata
        h5ad_paths: Iterable of paths to raw RNA h5ad files.
        output_dir: Path to folder to save outputs. It will have two sub-folders:
            "rna_qc_reports" will contain QC CSVs
            "pseudobulks" will contain pseudobulked h5ad and TSVs
        rna_qc_reports_dir: Path to folder to save RNA qc reports
        temp_pseudobulk_dir: Path to folder to write the per-pseudobulk AnnData pickles into
    Returns:
        defaultdict with keys being pseudobulk IDs, and values being a list of paths to pickled
            AnnData objects with raw RNA data corresponding to that pseudobulk ID
    """
    # write map of cell name to annotation
    _save_cell_name_to_annotation_mapping(metadata_loc=metadata_loc, output_dir=output_dir)
    # Iterate through h5ads
    futures = [
        executor.submit(
            _load_and_qc_h5ad,
            h5ad_path=h5ad_path,
            rna_qc_reports_dir=rna_qc_reports_dir,
            temp_pseudobulk_dir=temp_pseudobulk_dir,
        )
        for h5ad_path in h5ad_paths
    ]
    pseudobulk_adata_paths: defaultdict[PseudobulkName, list[Path]] = defaultdict(list)
    for adatas_result in as_completed(futures):
        for pseudobulk_id, adata_path in adatas_result.result():
            pseudobulk_adata_paths[pseudobulk_id].append(adata_path)

    return pseudobulk_adata_paths


def _aggregate_pseudobulk(
    pseudobulk_id: PseudobulkName,
    adata_paths: list[Path],
    rna_qc_reports_dir: Path,
    pseudobulked_rna_dir: Path,
) -> None:
    """Aggregate AnnDatas for this pseudobulk ID, save pseudobulked data, counts, and QC.

    Args:
        pseudobulk_id: ID for this pseudobulk
        adata_paths: paths to the pickled AnnData objects that correspond to this pseudobulk, as
            written by _load_and_qc_h5ad. They are removed once they have been read back in.
        rna_qc_reports_dir: Path to folder to save pseudobulked RNA QC reports
        pseudobulked_rna_dir: Path to folder to save pseudobulked RNA
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    gene_ref = _worker_data("gene_ref")
    with _worker_log_lock():
        logger.info(f"Aggregating pseudobulk_id: {pseudobulk_id}")
    adatas: list[ad.AnnData] = []
    for adata_path in adata_paths:
        adatas.append(ad.read_h5ad(f"{adata_path}", backed=False))
        adatas[-1].obs["pseudobulk_id"] = pseudobulk_id
        # the AnnData is in memory now, so free the temp folder as we go
        adata_path.unlink(missing_ok=True)
    p_qc = pd.concat([pd.DataFrame(adata.obs) for adata in adatas], axis=0)
    p_concat: ad.AnnData = ad.concat(adatas, axis=0)
    p_concat.var["gene_symbol"] = p_concat.var.index.map(gene_ref["gene_name"])

    # Save QC
    rna_qc_reports_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = rna_qc_reports_dir / f"{pseudobulk_id}.pseudobulked_cell_QC_metrics.tsv"
    p_qc.to_csv(f"{out_tsv}", sep="\t", index=False)
    # Save h5ad
    pseudobulked_rna_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad: Path = pseudobulked_rna_dir / f"{pseudobulk_id}.rna_counts_mtx.h5ad"
    p_concat.write(filename=f"{out_h5ad}")  # ty:ignore[missing-argument]
    # make pseudobulk
    counts_df_p = pd.DataFrame(p_concat.var.copy())
    counts_df_p["mt"] = counts_df_p.index.map(gene_ref["mt"]).fillna(False).astype(bool)
    counts_df_p["ribo"] = counts_df_p.index.map(gene_ref["ribo"]).fillna(False).astype(bool)
    counts_df_p["counts"] = p_concat.X.sum(axis=0).A1  # ty:ignore[unresolved-attribute]
    counts_df_p["CPM"] = (counts_df_p["counts"] / counts_df_p["counts"].sum()) * 1e6
    counts_df_p["log10CPM"] = np.log10(counts_df_p["CPM"] + 1)
    counts_df_p.to_csv(
        f"{pseudobulked_rna_dir}/{pseudobulk_id}.pseudobulk_expression.tsv.gz", sep="\t"
    )
    # as in _load_and_qc_h5ad: these AnnDatas are only reachable through reference cycles
    del adatas, p_concat
    gc.collect()


def _pseudobulk_rna_in_pool(
    *,
    executor: ProcessPoolExecutor,
    num_workers: int,
    metadata_loc: Path,
    input: Path,
    output_dir: Path,
    rna_qc_reports_dir: Path,
    temp_dir: Path,
    logger: logging.Logger,
) -> None:
    """Load, pseudobulk, and save the RNA data using the supplied pool of worker processes."""
    with tempfile.TemporaryDirectory(dir=f"{temp_dir}") as temp_pseudobulk_dir:
        logger.info(f"Loading and QC-ing h5ads with {num_workers} workers.")
        # Load raw RNA h5ads, save per-accession QC info, and group by pseudobulk ID
        pseudobulk_adata_paths = _load_and_qc_h5ads(
            executor=executor,
            metadata_loc=metadata_loc,
            h5ad_paths=_get_h5ad_paths(input, logger=logger),
            output_dir=output_dir,
            rna_qc_reports_dir=rna_qc_reports_dir,
            temp_pseudobulk_dir=Path(temp_pseudobulk_dir),
        )

        # aggregate across pseudobulks and save
        logger.info(
            f"Aggregating {len(pseudobulk_adata_paths)} pseudobulks with {num_workers} workers."
        )
        pseudobulked_rna_dir = output_dir / "pseudobulks"
        futures = [
            executor.submit(
                _aggregate_pseudobulk,
                pseudobulk_id=pseudobulk_id,
                adata_paths=adata_paths,
                rna_qc_reports_dir=rna_qc_reports_dir,
                pseudobulked_rna_dir=pseudobulked_rna_dir,
            )
            for pseudobulk_id, adata_paths in pseudobulk_adata_paths.items()
        ]
        for future in as_completed(futures):
            future.result()  # raise any exceptions in worker processes


def pseudobulk_rna(
    *,
    input: Path,
    output_dir: Path,
    metadata_loc: Path,
    gene_info: Path,
    num_workers: int = -1,
    temp_dir: Path = Path("."),
) -> None:
    """Separate RNA h5ad files by pseudorep.

    Args:
        input: Either a folder, in which case raw input RNA files will be globbed from this folder,
            or a file, with one input RNA file per line.
        output_dir: Path to folder to save outputs. It will have two sub-folders:
            "rna_qc_reports" will contain QC CSVs,
            "pseudobulks" will contain pseudobulked h5ad and TSVs.
        metadata_loc: Input annotations metadata file path.
        at_annotation_level: Whether to split pseudobulks only at the annotation level (if True) or
            also at the subsample level (if False).
        gene_info: Path to species-specific CSV of gene info.
        num_workers: Number of parallel workers to use. If <=0, use all available cores.
        temp_dir: Root folder for temporary files.
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    num_workers = num_workers if num_workers > 0 else cpu_count()
    # Load gene information and metadata once, and hand them to each worker process at start up
    gene_ref = utils.read_csv(gene_info, index_col=0)
    rna_qc_reports_dir = output_dir / "rna_qc_reports"
    # note the OOM kill count now, so that a worker that is killed later can be identified as having
    # run out of memory
    oom_kills_before = utils.oom_kill_count()
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_init_worker,
        initargs=(gene_ref, metadata_loc, Lock(), logger.getEffectiveLevel()),
    ) as executor:
        try:
            _pseudobulk_rna_in_pool(
                executor=executor,
                num_workers=num_workers,
                metadata_loc=metadata_loc,
                input=input,
                output_dir=output_dir,
                rna_qc_reports_dir=rna_qc_reports_dir,
                temp_dir=temp_dir,
                logger=logger,
            )
        except BrokenProcessPool:
            utils.exit_if_oom_killed(executor, oom_kills_before, logger)
            raise
