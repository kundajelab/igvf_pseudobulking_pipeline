import gc
import logging
import tempfile
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    as_completed,
    wait,
)
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import Lock, cpu_count
from pathlib import Path
from typing import Any, Final, cast

import anndata as ad
import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.sparse

from pseudobulk import sc_qc, utils
from pseudobulk.types import PseudobulkName

_LOGGER_NAME: Final[str] = "pseudobulk-rna"

# Reference data that every task needs, held once per worker process. Passing these as task
# arguments instead would re-pickle them (the metadata is ~300 MB in memory) for every h5ad.
# NOTE: pandas is not safe to use from several threads with the GIL disabled, so the workers are
# processes, and AnnData objects are handed between them as files rather than as arguments.
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


def _get_h5ad_paths(rna_input: Path, logger: logging.Logger) -> Iterable[Path]:
    """Get paths to raw RNA h5ad files from input.

    Args:
        rna_input: Either a folder, in which case raw input RNA files will be globbed from this
            folder, or a file, with one input RNA file per line.
        logger: Logger to output progress.

    Returns:
        Iterable of raw RNA h5ad files.
    """
    if rna_input.is_dir():
        logger.info(f"Finding h5ad files in {rna_input}")
        yield from rna_input.glob("*.h5ad")
    elif rna_input.is_file():
        logger.info(f"Reading input paths from file-of-files: {rna_input}")
        with rna_input.open("rt") as f_in:
            for h5ad_file in f_in:
                yield Path(h5ad_file.strip())
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
    # NOTE: equivalent to scanpy.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"],
    # percent_top=None, log1p=False, inplace=True), but works with backed AnnData
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
        obs.to_csv(
            f"{qc_report_path}",
            sep=sep,
            index=False,
        )
    return obs


_SPARSE_TYPES: Final = (
    scipy.sparse.csr_matrix,
    scipy.sparse.csc_matrix,
    scipy.sparse.csr_array,
    scipy.sparse.csc_array,
)
"""Types of sparse matrix that anndata may read from an h5ad."""

type _Rows = (
    scipy.sparse.csr_matrix
    | scipy.sparse.csc_matrix
    | scipy.sparse.csr_array
    | scipy.sparse.csc_array
    | np.ndarray
    | pd.DataFrame
)
"""Types of the AnnData elements that _read_rows can read rows of."""


def _read_rows(elem: h5py.Group | h5py.Dataset, row_idx: npt.NDArray[np.intp]) -> _Rows:
    """Read only the wanted rows of an AnnData element stored in an h5ad."""
    match elem.attrs.get("encoding-type"):
        # sparse matrices are stored as a group of their data, indices and indptr datasets
        case "csr_matrix" | "csc_matrix" if isinstance(elem, h5py.Group):
            rows = ad.io.sparse_dataset(elem)[row_idx]
            # NOTE: check for the sparse types rather than excluding the scalar that indexing a
            # single element gives: anndata annotates its sparse types in a way that type checkers
            # cannot resolve
            if not isinstance(rows, _SPARSE_TYPES):
                raise TypeError(f"selecting rows of {elem.name} in h5ad gave a {type(rows)}")
            return rows
        case "array" if isinstance(elem, h5py.Dataset):
            # h5py only fancy-indexes with a non-empty selection
            return elem[row_idx] if len(row_idx) > 0 else elem[:0]
        case _:
            # other encodings (e.g. dataframes in obsm) are rare and small: read them whole, then
            # subset
            match ad.io.read_elem(elem):
                case pd.DataFrame() as whole:
                    return whole.iloc[row_idx]
                case np.ndarray() as whole:
                    return np.take(whole, row_idx, axis=0)
                case _ as whole:
                    raise TypeError(
                        f"cannot select rows of {elem.name} in h5ad: unsupported type {type(whole)}"
                    )


def _save_adata_pseudobulk(
    h5ad_path: Path,
    pseudobulk_id: PseudobulkName,
    obs_path: Path,
    adata_path: Path,
) -> tuple[PseudobulkName, Path]:
    """Separate out the AnnData for the pseudobulk from a raw RNA h5ad, and save to adata_path.

    This is a separate task from _load_and_qc_h5ad, so that the pseudobulks of one large h5ad can be
    separated out in parallel. It reads the h5ad with h5py rather than anndata, because a backed
    AnnData still reads every layer whole into memory, and several of these tasks may be reading
    the same large h5ad at once.

    Args:
        h5ad_path: path to the raw RNA h5ad.
        pseudobulk_id: ID of the pseudobulk.
        obs_path: path to the pickled (row indices, obs) of the pseudobulk's cells, as written by
            _load_and_qc_h5ad. It is removed once it has been read back in.
        adata_path: path to save the pseudobulk's AnnData to.

    Returns:
        (pseudobulk ID, adata_path)
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    with _worker_log_lock():
        logger.info(f"\tIn {h5ad_path}: separating out AnnData for pseudobulk_id {pseudobulk_id}")
    row_idx, obs = cast(tuple[npt.NDArray[np.intp], pd.DataFrame], pd.read_pickle(obs_path))
    obs_path.unlink()
    with h5py.File(h5ad_path, "r") as h5ad:
        # NOTE: AnnData annotates obsm as holding only arrays and sequences, but like the layers it
        # also accepts sparse matrices and dataframes
        obsm: dict[str, Any] = {
            name: _read_rows(elem, row_idx)
            for name, elem in cast(h5py.Group, h5ad.get("obsm", {})).items()
        }
        # obsp and varm are not read: ad.concat drops them when the pseudobulk is aggregated
        pseudobulk_adata = ad.AnnData(
            X=_read_rows(h5ad["X"], row_idx),
            obs=obs,
            var=cast(pd.DataFrame, ad.io.read_elem(h5ad["var"])),
            layers={
                name: _read_rows(elem, row_idx)
                for name, elem in cast(h5py.Group, h5ad.get("layers", {})).items()
            },
            obsm=obsm,
        )
    ad.io.write_h5ad(adata_path, pseudobulk_adata)
    # garbage collection is needed because of how AnnData objects are stored.
    del pseudobulk_adata
    gc.collect()
    return pseudobulk_id, adata_path


def _load_and_qc_h5ad(
    h5ad_path: Path,
    rna_qc_reports_dir: Path,
    temp_pseudobulk_dir: Path,
) -> list[tuple[PseudobulkName, Path, Path]]:
    """Load one raw RNA h5ad, save its QC, and write out the obs of each pseudobulk ID's cells.

    Returns:
        list of (pseudobulk ID, path to the pickled (row indices, obs) for that pseudobulk, path to
            save the pseudobulk's AnnData to), to pass to _save_adata_pseudobulk.
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
    obs = pd.DataFrame(adata.obs)
    adata.file.close()
    # AnnData objects take part in reference cycles, so nothing but the cycle collector can free
    # them, and it is triggered by the number of allocations rather than by their size. A handful
    # of hundred-MB objects never reaches that threshold, so without collecting here each task
    # leaves its whole h5ad behind, growing with each file without bound, until the task was killed
    # for running out of memory.
    del adata
    gc.collect()

    # Save the obs of each pseudobulk ID's cells
    pseudobulk_paths: list[tuple[PseudobulkName, Path, Path]] = []
    # observed=True because pseudobulk_id is a categorical whose categories span every
    # pseudobulk in the run: the default would yield a group for each of them for every
    # accession, almost all empty, and write an h5ad for each one.
    for pseudobulk_id, pseudobulk_metadata in metadata_df.groupby(
        "pseudobulk_id", sort=False, group_keys=False, as_index=False, observed=True
    ):
        row_idx = np.flatnonzero(obs.index.isin(set(pseudobulk_metadata["barcode_sample"])))
        obs_path = temp_pseudobulk_dir / f"{pseudobulk_id}.{accession}.obs.pkl"
        pd.to_pickle((row_idx, obs.iloc[row_idx, :]), obs_path)
        pseudobulk_paths.append(
            (
                PseudobulkName(f"{pseudobulk_id}"),
                obs_path,
                temp_pseudobulk_dir / f"{pseudobulk_id}.{accession}.h5ad",
            )
        )
    return pseudobulk_paths


def _load_and_qc_h5ads(
    executor: ProcessPoolExecutor,
    h5ad_paths: Iterable[Path],
    rna_qc_reports_dir: Path,
    temp_pseudobulk_dir: Path,
) -> defaultdict[PseudobulkName, list[Path]]:
    """Load raw RNA h5ads and separate into pseudobulked AnnData objects.

    Also save observation and QC data for each accession ID into rna_qc_reports_dir.

    Args:
        executor: process pool to load the h5ads with. Its workers must have been initialized with
            _init_worker, which holds the metadata and gene reference data they need.
        h5ad_paths: Iterable of paths to raw RNA h5ad files.
        rna_qc_reports_dir: Path to folder to save RNA qc reports
        temp_pseudobulk_dir: Path to folder to write the temporary per-pseudobulk files into

    Returns:
        defaultdict with keys being pseudobulk IDs, and values being a list of paths to temporary
            AnnData objects with raw RNA data corresponding to that pseudobulk ID
    """
    # QC each h5ad, then as each finishes, separate out its pseudobulks as their own tasks, so that
    # a single large h5ad is still split up across the workers
    qc_futures: dict[Future[list[tuple[PseudobulkName, Path, Path]]], Path] = {
        executor.submit(
            _load_and_qc_h5ad,
            h5ad_path=h5ad_path,
            rna_qc_reports_dir=rna_qc_reports_dir,
            temp_pseudobulk_dir=temp_pseudobulk_dir,
        ): h5ad_path
        for h5ad_path in h5ad_paths
    }
    save_futures: set[Future[tuple[PseudobulkName, Path]]] = set()
    pseudobulk_adata_paths: defaultdict[PseudobulkName, list[Path]] = defaultdict(list)
    while qc_futures or save_futures:
        done, _ = wait([*qc_futures, *save_futures], return_when=FIRST_COMPLETED)
        for future in done:
            h5ad_path = qc_futures.pop(
                cast(Future[list[tuple[PseudobulkName, Path, Path]]], future), None
            )
            if h5ad_path is None:
                save_future = cast(Future[tuple[PseudobulkName, Path]], future)
                save_futures.remove(save_future)
                pseudobulk_id, adata_path = save_future.result()
                pseudobulk_adata_paths[pseudobulk_id].append(adata_path)
            else:
                save_futures.update(
                    executor.submit(
                        _save_adata_pseudobulk,
                        h5ad_path=h5ad_path,
                        pseudobulk_id=pseudobulk_id,
                        obs_path=obs_path,
                        adata_path=adata_path,
                    )
                    for pseudobulk_id, obs_path, adata_path in future.result()
                )

    return pseudobulk_adata_paths


def _aggregate_pseudobulk(
    pseudobulk_id: PseudobulkName,
    adata_paths: list[Path],
    rna_qc_reports_dir: Path,
    pseudobulked_rna_dir: Path,
) -> int:
    """Aggregate AnnDatas for this pseudobulk ID, save pseudobulked data, counts, and QC.

    Args:
        pseudobulk_id: ID for this pseudobulk
        adata_paths: paths to the temporary AnnData h5ads that correspond to this pseudobulk, as
            written by _save_adata_pseudobulk. They are removed once they have been read back in.
        rna_qc_reports_dir: Path to folder to save pseudobulked RNA QC reports
        pseudobulked_rna_dir: Path to folder to save pseudobulked RNA
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    gene_ref = _worker_data("gene_ref")
    with _worker_log_lock():
        logger.info(f"Aggregating pseudobulk_id: {pseudobulk_id}")
    adatas: list[ad.AnnData] = []
    # adata_paths are in the order their h5ads finished, so sort them (by accession, as they are all
    # for this pseudobulk) to keep the order of cells, and so the outputs, reproducible
    for adata_path in sorted(adata_paths):
        adatas.append(ad.read_h5ad(f"{adata_path}", backed=False))
        adatas[-1].obs["pseudobulk_id"] = pseudobulk_id
        # the AnnData is in memory now, so free the temp folder as we go
        adata_path.unlink(missing_ok=True)
    p_qc = pd.concat([pd.DataFrame(adata.obs) for adata in adatas], axis=0)
    p_concat: ad.AnnData = ad.concat(adatas, axis=0)
    p_concat.var["gene_symbol"] = p_concat.var.index.map(gene_ref["gene_name"])
    num_pseudobulk_rows = len(p_concat)
    if num_pseudobulk_rows == 0:
        with _worker_log_lock():
            logger.warning(f"pseudobulk {pseudobulk_id} is empty.")

    # Save QC
    rna_qc_reports_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = rna_qc_reports_dir / f"{pseudobulk_id}.pseudobulked_cell_QC_metrics.tsv"
    p_qc.to_csv(f"{out_tsv}", sep="\t", index=False)
    # Save h5ad
    pseudobulked_rna_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad: Path = pseudobulked_rna_dir / f"{pseudobulk_id}.rna_counts_mtx.h5ad"
    p_concat.write(filename=f"{out_h5ad}")
    # make pseudobulk
    counts_df_p = pd.DataFrame(p_concat.var.copy())
    counts_df_p["mt"] = counts_df_p.index.map(gene_ref["mt"]).fillna(False).astype(bool)
    counts_df_p["ribo"] = counts_df_p.index.map(gene_ref["ribo"]).fillna(False).astype(bool)
    counts_df_p["counts"] = p_concat.X.sum(axis=0, dtype=np.uint64).A1.astype(np.uint64)  # ty:ignore[unresolved-attribute]
    counts_df_p["CPM"] = (counts_df_p["counts"] / counts_df_p["counts"].sum()) * 1.0e6
    # accurately calculate log10(1 + CPM), then round to 14 digits of accuracy to keep repeatable
    counts_df_p["log10CPM"] = np.around(np.log1p(counts_df_p["CPM"]) / np.log(10.0), 14)
    counts_df_p.to_csv(
        f"{pseudobulked_rna_dir}/{pseudobulk_id}.pseudobulk_expression.tsv.gz",
        sep="\t",
        compression=utils.COMPRESSION_DICT,
    )
    # as in _load_and_qc_h5ad: these AnnDatas are only reachable through reference cycles
    del adatas, p_concat
    gc.collect()
    return num_pseudobulk_rows


def _pseudobulk_rna_in_pool(
    *,
    executor: ProcessPoolExecutor,
    num_workers: int,
    rna_input: Path,
    output_dir: Path,
    rna_qc_reports_dir: Path,
    temp_dir: Path,
    logger: logging.Logger,
) -> int:
    """Load, pseudobulk, and save the RNA data using the supplied pool of worker processes."""
    with tempfile.TemporaryDirectory(dir=f"{temp_dir}") as temp_pseudobulk_dir:
        logger.info(f"Loading and QC-ing h5ads with {num_workers} workers.")
        # Load raw RNA h5ads, save per-accession QC info, and group by pseudobulk ID
        pseudobulk_adata_paths = _load_and_qc_h5ads(
            executor=executor,
            h5ad_paths=_get_h5ad_paths(rna_input, logger=logger),
            rna_qc_reports_dir=rna_qc_reports_dir,
            temp_pseudobulk_dir=Path(temp_pseudobulk_dir),
        )

        # aggregate across pseudobulks and save
        logger.info(
            f"Aggregating {len(pseudobulk_adata_paths)} pseudobulks with {num_workers} workers."
        )
        pseudobulked_rna_dir = output_dir / "pseudobulks"
        max_pseudobulk_rows = 0
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
            # raise any exceptions in worker processes, and collect number of pseudobulk rows
            max_pseudobulk_rows = max(max_pseudobulk_rows, future.result())
        return max_pseudobulk_rows


def pseudobulk_rna(
    *,
    rna_input: Path,
    output_dir: Path,
    metadata_loc: Path,
    gene_info: Path,
    num_workers: int = -1,
    temp_dir: Path = Path("."),
) -> None:
    """Separate RNA h5ad files by pseudobulk.

    Args:
        rna_input: Either a folder, in which case raw input RNA files will be globbed from this
            folder, or a file, with one input RNA file per line.
        output_dir: Path to folder to save outputs. It will have two sub-folders:
            "rna_qc_reports" will contain QC CSVs,
            "pseudobulks" will contain pseudobulked h5ad and TSVs.
        metadata_loc: Input annotations metadata file path.
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
            num_pseudobulk_rows = _pseudobulk_rna_in_pool(
                executor=executor,
                num_workers=num_workers,
                rna_input=rna_input,
                output_dir=output_dir,
                rna_qc_reports_dir=rna_qc_reports_dir,
                temp_dir=temp_dir,
                logger=logger,
            )
        except BrokenProcessPool:
            utils.exit_if_oom_killed(executor, oom_kills_before, logger)
            raise

    if num_pseudobulk_rows == 0:
        raise RuntimeError("No non-empty pseudobulks were found.")
