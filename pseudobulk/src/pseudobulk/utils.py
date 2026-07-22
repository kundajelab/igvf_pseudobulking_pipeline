import re
import signal
import unicodedata
from collections import defaultdict
from collections.abc import (
    Hashable,
    Iterable,
    Mapping,
    Sequence,
)
from concurrent.futures import ProcessPoolExecutor
from contextlib import AbstractContextManager, contextmanager, nullcontext
from logging import Logger
from pathlib import Path
from types import MappingProxyType
from typing import (
    Callable,
    Final,
    Generator,
    Literal,
    TextIO,
    TypeAlias,
    overload,
)

import numpy as np
import pandas as pd
import scipy.sparse

from pseudobulk.barcode_qc import BarcodeQc
from pseudobulk.tss import update_insertions_range
from pseudobulk.types import (
    COUNTS_ARRAY,
    COUNTS_MATRIX,
    POS_ARRAY,
    POS_DTYPE,
    Barcode,
    Contig,
    PseudobulkName,
)

LogLock: TypeAlias = AbstractContextManager[bool]
"""A lock used to keep log messages from interleaving.

Tools that parallelize with processes rather than threads must supply a multiprocessing lock, since
a threading lock only serializes the threads within one process.
"""

OOM_EXIT_STATUS: Final[int] = 137
"""Exit status of a process killed for running out of memory (128 + SIGKILL)."""

_CGROUP_ROOT: Final[Path] = Path("/sys/fs/cgroup")
"""Mount point of the cgroup filesystem."""

_CGROUP_OOM_FILE_NAMES: Final[tuple[str, ...]] = (
    "memory.events",  # cgroup v2
    "memory.oom_control",  # cgroup v1, but only reports oom_kill on kernel 4.13 and later
)
"""Files reporting the OOM kill count of a cgroup, in the order they should be tried."""

ATAC_QC_COLUMNS: Final[tuple[str, ...]] = tuple(BarcodeQc.header_columns())
"""Expected columns in ATAC QC files."""

RNA_QC_COLUMNS: Final[tuple[str, ...]] = (
    "analysis_set_accession",
    "barcode_sample",
    "annotated",
    "pseudobulk_id",
    "rna_read_count",
    "gene_count",
    "pct_mito",
    "pct_ribo",
)
"""Expected columns in RNA QC files."""

FRONT_QC_COLUMNS: Final[tuple[str, ...]] = (
    "analysis_set_accession",
    "barcode_sample",
    "annotated",
    "found_in_rna",
    "found_in_atac",
    "pseudobulk_id",
)
"""Columns that should come first in combined RNA + ATAC QC."""

MANDATORY_METADATA_COLS: Final[tuple[str, ...]] = (
    "barcode_sample",
    "cell_name",
    "cell_description",
    "CL_id",
    "CL_term_name",
    "subsample",
    "analysis_set_accession",
)
"""Columns that must be present in annotations files."""

OPTIONAL_METADATA_COLS: Final[tuple[str, ...]] = (
    "matrix_file_accession",
    "fragments_file_accession",
    "cell_qualifier",
)
"""Columns that may be present in annotations files and should be used if they are."""

ADDED_METADATA_COLS: Final[dict[str, tuple[str, ...]]] = {
    "annotation": ("cell_name",),
    "pseudobulk_id": ("cell_name", "subsample"),
}
"""Usedful derived columns that are not in the raw annotations file, and their raw dependencies."""

METADATA_COLS: Final[tuple[str, ...]] = (
    MANDATORY_METADATA_COLS + OPTIONAL_METADATA_COLS + tuple(ADDED_METADATA_COLS.keys())
)
"""All potentially wanted columns"""


def _own_cgroup_dirs() -> list[Path]:
    """Get this process's own cgroup directories, closest first, then their ancestors.

    The OOM kill may be reported against the cgroup of the process or against any of its ancestors
    (under SLURM the limit is set on the job's cgroup, not on the root one).
    """
    cgroup_dirs: list[Path] = []
    try:
        cgroup_report = Path("/proc/self/cgroup").read_text()
    except OSError:
        return [_CGROUP_ROOT]
    for line in cgroup_report.splitlines():
        # cgroup v2 lines are "0::<path>", cgroup v1 lines are "<id>:<controllers>:<path>"
        _hierarchy, _, rest = line.partition(":")
        controllers, _, cgroup_path = rest.partition(":")
        if controllers and "memory" not in controllers.split(","):
            continue
        base = _CGROUP_ROOT if not controllers else _CGROUP_ROOT / "memory"
        cgroup_dir = base / cgroup_path.lstrip("/")
        while True:
            cgroup_dirs.append(cgroup_dir)
            if cgroup_dir == base or base not in cgroup_dir.parents:
                break
            cgroup_dir = cgroup_dir.parent
    cgroup_dirs.append(_CGROUP_ROOT)
    return cgroup_dirs


def oom_kill_count() -> int | None:
    """Get the number of times the kernel has OOM-killed a process in this cgroup.

    Returns:
        The OOM kill count, or None if it cannot be determined. That is the case on macOS (no
        cgroups at all) and on cgroup v1 kernels before 4.13, such as CentOS 7, which report
        whether a cgroup is currently out of memory but not how many processes have been killed.
    """
    for cgroup_dir in _own_cgroup_dirs():
        for oom_file_name in _CGROUP_OOM_FILE_NAMES:
            try:
                oom_report = (cgroup_dir / oom_file_name).read_text()
            except OSError:
                continue
            for line in oom_report.splitlines():
                key, _, value = line.partition(" ")
                if key == "oom_kill":
                    try:
                        return int(value)
                    except ValueError:
                        return None
    return None


def killed_worker_count(executor: ProcessPoolExecutor, logger: Logger) -> int | None:
    """Get how many of a broken pool's worker processes were killed with SIGKILL.

    SIGKILL is what the kernel OOM killer sends, and nothing else in normal operation sends it to a
    worker: a worker that crashes exits on its own signal (e.g. SIGSEGV) or status, and the pool
    shuts its remaining workers down with SIGTERM. So a SIGKILLed worker means the process was
    killed from outside, which on a compute node means it ran out of memory.

    Args:
        executor: the process pool whose workers died.
    Returns:
        The number of workers killed with SIGKILL, or None if the pool does not expose its worker
        processes (it is a private attribute, so treat it as unavailable rather than assuming zero).
    """
    worker_processes = getattr(executor, "_processes", None)
    if worker_processes is None:
        logger.warning("Unable to get worker processes form executor._processes")
        return None
    exit_codes = defaultdict(lambda: 0)
    for process in worker_processes.values():
        exit_codes[process.exitcode] += 1
    exit_codes_str = ", ".join(f"{code}: {counts}" for code, counts in exit_codes.items())
    logger.warning(f"Counts of exit codes from workers: {exit_codes_str}")
    shutdown_codes = {-signal.SIGKILL, -signal.SIGTERM}
    return sum(
        1 if process.exitcode in shutdown_codes else 0 for process in worker_processes.values()
    )


def exit_if_oom_killed(
    executor: ProcessPoolExecutor, oom_kills_before: int | None, logger: Logger
) -> None:
    """Exit with the out-of-memory status if a worker process was killed for running out of memory.

    A worker process that is OOM-killed only reaches the parent process as a BrokenProcessPool, so
    without this the tool would exit with a generic failure status and the pipeline could not tell
    that the task should be retried with more memory. Only exits when the kill can actually be
    attributed to running out of memory: otherwise it returns, and the caller should re-raise the
    original error so that a crashing worker is not misreported (and retried) as an OOM.

    Args:
        executor: the process pool whose workers died.
        oom_kills_before: OOM kill count from oom_kill_count(), taken before the work started.
        logger: Logger to report the out-of-memory kill with.
    """
    match killed_worker_count(executor, logger):
        case 0:
            # got a definitive count, but no workers were killed
            return
        case None:
            # did not get a difinitive count, try to get a reporting from cgroup counter
            # for linux kernels new enough to report it
            oom_kills_after = oom_kill_count()
            if oom_kills_before is None or oom_kills_after is None:
                logger.error(
                    "Unable to determine if any workers were killed due to missing cgroup data."
                )
                return
            elif oom_kills_after > oom_kills_before:
                logger.error(
                    f"the kernel reported {oom_kills_after - oom_kills_before} out-of-memory kill(s) "
                    "in this cgroup"
                )
                raise SystemExit(OOM_EXIT_STATUS)
            else:
                # got a definitive count, but no oom kills:
                return
        case _ as num_killed:
            logger.error(
                f"{num_killed} worker processes ran out of memory. Exiting with status {OOM_EXIT_STATUS} so "
                "that the task is retried with more memory."
            )
            raise SystemExit(OOM_EXIT_STATUS)


def load_metadata(
    metadata_loc: Path,
    wanted_cols: Iterable[str] | None = METADATA_COLS,
) -> pd.DataFrame:
    """Load metadata file and perform basic checks."""
    # determine which columns we need to read in
    if wanted_cols is None:
        usecols = None
        added_cols = set()
    else:
        usecols = []
        dependency_cols = set()
        added_cols = set()
        available_columns = (
            read_csv(metadata_loc, nrows=0).columns
            if len(set(wanted_cols).intersection(OPTIONAL_METADATA_COLS)) > 0
            else None
        )

        for col in wanted_cols:
            col_dependencies = ADDED_METADATA_COLS.get(col, None)
            if col_dependencies is None:
                if available_columns is None or col in available_columns:
                    usecols.append(col)
            else:
                added_cols.add(col)
                dependency_cols.update(col_dependencies)
        usecols.extend(dependency_cols.difference(usecols))
    # read in the TSV
    df = read_csv(metadata_loc, usecols=usecols)
    if len(added_cols) > 0:
        # if we need to add any derived cols, add them
        _add_derived_metadata(df, added_cols=added_cols)
    if wanted_cols is not None:
        # drop the source columns that were only read so that the added ones could be derived
        df = df.loc[:, [col for col in wanted_cols if col in df.columns]]
    # check 'subsample' contains valid values
    if "subsample" in df.columns:
        for subsample in df["subsample"].unique():
            if "-" in subsample:
                raise ValueError(
                    f"'subsample' column contained '{subsample}' containing invalid annotation values"
                )
    # save space by converting repeat values to categoricals
    for col_name, col_values in df.items():
        if col_values.nunique() * 2 < len(col_values):
            df[col_name] = col_values.astype("category")
    return df


def sanitize_to_ascii_underscore(text: str) -> str:
    """Replace unicode characters with similar ascii, replace whitespace and commas with underscores."""
    # 1. Normalize Unicode to NFKD form to separate characters from accents
    # 2. Encode to ASCII and ignore characters that cannot be converted
    # 3. Decode back to a string
    text = unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")

    # 4. Replace one or more commas or whitespace characters with a single underscore
    # \s+ matches spaces, tabs, and newlines
    return re.sub(r"[,\s]+", "_", text)


def _get_map_to_sanitized(cell_name: pd.Series) -> dict[str, str]:
    """Get a 1-to-1 map from cell_name to safe ascii annotation.

    Append indices as needed to guarantee uniqueness.
    """
    sorted_cell_names = sorted(cell_name.unique())
    reversed_dict = defaultdict(list)
    for cell_name in sorted_cell_names:
        reversed_dict[sanitize_to_ascii_underscore(cell_name)].append(cell_name)
    while any(len(cell_names) > 1 for cell_names in reversed_dict.values()):
        old_reversed_dict, reversed_dict = reversed_dict, defaultdict(list)
        for sanitized_name in sorted(old_reversed_dict.keys()):
            cell_names = old_reversed_dict[sanitized_name]
            if len(cell_names) == 1:
                reversed_dict[sanitized_name].append(cell_names[0])
            else:
                for idx, cell_name in enumerate(sorted(cell_names), start=1):
                    reversed_dict[f"{sanitized_name}_{idx}"].append(cell_name)
    return {cell_names[0]: sanitized_name for sanitized_name, cell_names in reversed_dict.items()}


def _add_derived_metadata(metadata_df, added_cols: set[str]) -> None:
    """
    Update metadata_df in place with new 'annotation' and 'pseudobulk_id' columns if they are in added_cols
    """
    # Cell name to annotation mapping
    cell_name_to_annotation_dict = _get_map_to_sanitized(metadata_df["cell_name"])
    # Map cell names to annotations in metadata_df
    metadata_df["annotation"] = (
        metadata_df["cell_name"].map(cell_name_to_annotation_dict).astype("category")
    )
    # subsample can have commas in it sometimes, so we need to sanitize that too
    subsample_dict = _get_map_to_sanitized(metadata_df["subsample"])
    if "pseudobulk_id" in added_cols:
        metadata_df["pseudobulk_id"] = metadata_df.apply(
            lambda row: f"{row['annotation']}-{subsample_dict[row['subsample']]}", axis=1
        ).astype("category")


def map_barcodes_to_pseudobulks(metadata_df: pd.DataFrame) -> Mapping[Barcode, PseudobulkName]:
    """Get mapping from barcodes to pseudobulk IDs."""
    return MappingProxyType(
        {
            Barcode(row["barcode_sample"]): PseudobulkName(row["pseudobulk_id"])
            for _, row in metadata_df.iterrows()
        }
    )


def load_tss_locs(tss: Path) -> dict[Contig, tuple[POS_ARRAY, POS_ARRAY]]:
    """Load Transcription Start Site (TSS) locations.

    Args:
        tss: Path to TSS locations file. The file should be a tab-separated values (TSV)
            file with columns "gene", "transcript", "chro", "TSS", and "strand".
    Returns:
        A dictionary mapping chromosome names to tuples of numpy arrays. Each tuple contains two
        numpy arrays: 1. TSS positions on that chromosome, 2. strand (+1 or -1)
    """
    tss_locs_df = read_csv(tss)
    tss_by_chr_np: dict[Contig, tuple[np.ndarray, np.ndarray]] = dict()
    for chro, chro_tss_df in tss_locs_df.groupby("chro", sort=False):
        tss_positions = chro_tss_df["TSS"].values - 1  # GTF is 1-based but fragments are 0-based
        strand_signs = np.where(chro_tss_df["strand"] == "+", 1, -1)
        df = (
            pd.DataFrame({"TSS": tss_positions, "strand": strand_signs})
            .sort_values(by=["TSS", "strand"])
            .drop_duplicates()
        )
        tss_by_chr_np[Contig(chro)] = (  # ty:ignore[invalid-argument-type]
            df["TSS"].to_numpy().astype(POS_DTYPE),
            df["strand"].to_numpy().astype(POS_DTYPE),
        )
    return tss_by_chr_np


def update_tss_insertions(
    tss_insertions: COUNTS_ARRAY,
    position: POS_DTYPE.type,
    tss_vec: POS_ARRAY,
    strand_vec: POS_ARRAY,
) -> None:
    """Update TSS counts with overlapping Transcription Start Site (TSS) in the given TSS vectors.

    Args:
        ts_insertions: Numpy array of counts of TSS.
        position: Genomic position to check (0-based).
        tss_vec: Numpy array of TSS positions on the same chromosome (0-based).
        strand_vec: Numpy array of strand signs corresponding to the TSS positions (+1 or -1).
    Returns:
        Numpy array of distances from each TSS to position (in the strand direction). If there are
        no TSS within half_window of the position, returns an empty array.
    """
    # find index where overlaps start on the left side
    half_window = POS_DTYPE.type(len(tss_insertions) // 2)
    tss_left_idx = np.searchsorted(tss_vec, position - half_window, side="left")
    # find relative (to left) index of end of overlaps on the right side
    tss_idx_delta = np.searchsorted(tss_vec[tss_left_idx:], position + half_window, side="right")
    if tss_idx_delta > 0:
        # there is >= 1 overlap, update the counts of insertion sites
        tss_right_idx = tss_left_idx + tss_idx_delta
        update_insertions_range(
            tss_insertions, tss_vec, strand_vec, tss_left_idx, tss_right_idx, half_window, position
        )


@contextmanager
def create_and_write(path: Path, mode: Literal["wt", "at"] = "wt") -> Generator[TextIO, None, None]:
    """Open path for writing, creating parent as necessary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, mode=mode) as f_out:
        yield f_out


def elapsed_time(elapsed_secs: float) -> str:
    """Convert float time in seconds to human-readable string."""
    minutes, seconds = divmod(elapsed_secs, 60.0)
    minutes = int(minutes)
    if minutes == 0:
        return f"{seconds:.3f}s"
    hours, minutes = divmod(minutes, 60)
    if hours == 0:
        return f"{minutes}:{seconds:06.3f}"
    else:
        return f"{hours}:{minutes:02d}:{seconds:06.3f}"


def get_sep_from_path(csv_or_tsv: Path) -> Literal["\t", ","]:
    """Guess the file separator from the path."""
    match csv_or_tsv.suffixes:
        case [*_parts, ".tsv"] | [*_parts, ".tsv", ".gz"]:
            return "\t"
        case [*_parts, ".csv"] | [*_parts, ".csv", ".gz"]:
            return ","
        case _:
            raise ValueError(
                f"Could not infer separator for file {csv_or_tsv} with suffix {csv_or_tsv.suffix}. Please "
                "specify separator with 'sep' argument."
            )


def read_csv(path: Path, sep: str | None = None, **kwargs) -> pd.DataFrame:
    """Read (possibly compressed) CSV or TSV using suffix to determine field separator."""
    if sep is None:
        sep = get_sep_from_path(path)
    return pd.read_csv(
        f"{path}",
        sep=sep,
        engine="c",
        low_memory=False,
        **kwargs,
    )


@overload
def load_atac_qc(
    atac_qc_dir: Path,
    need_tss: Literal[False],
    logger: Logger,
    log_lock: LogLock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
    accessions: set[str] | None = None,
) -> pd.DataFrame: ...


@overload
def load_atac_qc(
    atac_qc_dir: Path,
    need_tss: Literal[True],
    logger: Logger,
    log_lock: LogLock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
    accessions: set[str] | None = None,
) -> tuple[pd.DataFrame, COUNTS_MATRIX]: ...


def load_atac_qc(
    atac_qc_dir: Path,
    need_tss: bool,
    logger: Logger,
    log_lock: LogLock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
    accessions: set[str] | None = None,
) -> tuple[pd.DataFrame, COUNTS_MATRIX] | pd.DataFrame:
    """
    Load combined atac qc (generated per analysis accession, not by pseudobulk).

    If there are no ATAC QC files, return empty DataFrame and counts matrix.

    Args:
        atac_qc_dir: Path to folder with ATAC-seq QC files (should be TSVs and sparse matrices)
        need_tss: If True, load the Transcription Start Sites (TSS) sparse matrices.
        logger: Logger to output progress
        log_lock: An optional lock to use when logging. If None, no lock is used.
        row_filter: An optional func from an ATAC QC DataFrame to a series of bool describing if
            each row is wanted. If provided filter DataFrames (and TSS matrices) to those rows. If
            None, keep all rows.
        usecols: An optional selector of which columns to use
            if a sequence of strings, columns with those names are kept.
            if a sequence of indices, those column indices are kept.
            if a func from column label to bool, columns where the func return True are kept.
            If None, keep all columns.
        accessions: An optional set of analysis set accessions to load QC for. If provided, only the
            QC files of those accessions are read, which lets a worker load just the QC it needs. If
            None, load the QC of every accession in atac_qc_dir.
    Returns:
        combined ATAC-seq QC DataFrame and sparse TSS counts matrix if need_tss is True,
        otherwise just combined ATAC-seq DataFrame
    """
    _log_lock = log_lock if log_lock is not None else nullcontext()
    with _log_lock:
        logger.info(f"Loading ATAC QC from {atac_qc_dir}")
    atac_combined_qc_list: list[pd.DataFrame] = []
    tss_list: list[COUNTS_MATRIX] = []
    for tsv in atac_qc_dir.glob("*.tsv"):
        analysis_set_accession = tsv.name.split(".", 1)[0]
        if accessions is not None and analysis_set_accession not in accessions:
            continue
        with _log_lock:
            logger.info(f"loading QC for analysis accession {analysis_set_accession}")
        atac_qc = read_csv(tsv, usecols=usecols)
        if row_filter is None:
            row_is_wanted = None
        else:
            row_is_wanted = row_filter(atac_qc)
            atac_qc = atac_qc.loc[row_is_wanted]
        atac_combined_qc_list.append(atac_qc)
        matrix_file = atac_qc_dir / f"{analysis_set_accession}_tss_matrix.npz"
        if need_tss:
            tss = scipy.sparse.load_npz(f"{matrix_file}")
            if row_is_wanted is not None:
                tss = tss[row_is_wanted.to_numpy(), :]
            tss_list.append(tss)
    if len(atac_combined_qc_list) == 0:
        # there was no ATAC QC present, return empty objects
        atac_cols = ATAC_QC_COLUMNS
        match usecols:
            case func if callable(usecols):
                atac_cols = [col for col in atac_cols if func(col)]
            case [col, *_cols] if isinstance(col, str):
                atac_cols = list(usecols)
            case [col, *_cols] if isinstance(col, int):
                atac_cols = [atac_cols[idx] for idx in usecols]  # ty:ignore[invalid-argument-type]
            case _:
                pass

        atac_combined_qc = pd.DataFrame([], columns=pd.Index(atac_cols))
        tss_counts = scipy.sparse.csr_array(np.empty((0, 4001), np.uint16))
    else:
        atac_combined_qc = pd.concat(atac_combined_qc_list, axis=0, ignore_index=True)
        tss_counts = scipy.sparse.vstack(tss_list) if need_tss else None

    return (atac_combined_qc, tss_counts) if need_tss else atac_combined_qc


def _reorder_qc_columns(combined_qc: pd.DataFrame) -> pd.DataFrame:
    """Reorder columns to move selected columns to the front."""
    # Move shared columns to the front, and fill out the remaining columns in order afterwards
    col_order = [x for x in FRONT_QC_COLUMNS if x in combined_qc.columns] + [
        x for x in combined_qc.columns if x not in set(FRONT_QC_COLUMNS)
    ]
    # do the reorder
    return combined_qc[col_order]


def merge_rna_and_atac_qc(
    identifier: str,
    rna_qc: pd.DataFrame,
    atac_qc: pd.DataFrame,
    logger: Logger,
    log_lock: LogLock | None = None,
) -> pd.DataFrame:
    """Combined RNA QC and ATAC QC into one DataFrame object by merging on shared columns.

    Handles cases where DataFrames are empty, but requires columns to exist and be correct.

    Args:
        identifier: String to identify data being merged (e.g. accession set ID or pseudobulk ID).
        rna_qc: DataFrame
        logger: Logger to output progress
        log_lock: An optional lock to use when logging. If None, no lock is used
    Returns:
        Combined QC DataFrame with columns from both RNA and ATAC QC.
    """
    _log_lock = log_lock if log_lock is not None else nullcontext()
    if len(rna_qc.columns) == 0:
        raise ValueError("RNA QC column labels must be supplied, even if RNA QC is empty.")
    if len(atac_qc.columns) == 0:
        raise ValueError("ATAC QC column labels must be supplied, even if ATAC QC is empty.")
    match len(rna_qc), len(atac_qc):
        case 0, 0:
            # RNA and ATAC QC are empty (but contain the correct columns), return empty DataFrame
            combined_qc = pd.DataFrame([], columns=rna_qc.columns + atac_qc.columns)
            with _log_lock:
                logger.info(f"No RNA QC or ATAC QC for {identifier}")
        case 0, _:
            # RNA QC is empty (but will contain the correct columns)
            combined_qc = atac_qc.copy()
            for col in rna_qc.columns:
                if col in combined_qc.columns:
                    continue
                combined_qc[col] = float("NaN")
            with _log_lock:
                logger.info(f"No RNA QC for {identifier}")
        case _, 0:
            # ATAC QC is missing (but contains the correct columns)
            combined_qc = rna_qc.copy()
            for col in atac_qc.columns:
                if col in combined_qc.columns:
                    continue
                combined_qc[col] = float("NaN")
            with _log_lock:
                logger.info(f"No ATAC QC for {identifier}")
        case _:
            # both present, merge. Drop duplicate pseudobulk_id from atac_qc

            # Confirm that ATAC and RNA cells match
            combined_qc = pd.merge(
                rna_qc,
                atac_qc.drop(columns="pseudobulk_id"),
                how="outer",
                on=["analysis_set_accession", "barcode_sample", "annotated"],
            )
            with _log_lock:
                logger.info(f"Merged RNA QC and ATAC QC for {identifier}")

    # NOTE: fillna(0.) instead of fillna(False) to avoid a deprecation warning
    if "found_in_rna" in combined_qc.columns:
        combined_qc["found_in_rna"] = combined_qc["found_in_rna"].fillna(0.0).astype(bool)
    if "found_in_atac" in combined_qc.columns:
        combined_qc["found_in_atac"] = combined_qc["found_in_atac"].fillna(0.0).astype(bool)
    if "num_frags" in combined_qc.columns:
        combined_qc["num_frags"] = combined_qc["num_frags"].astype(pd.UInt64Dtype())
    if "rna_read_count" in combined_qc.columns:
        combined_qc["rna_read_count"] = combined_qc["rna_read_count"].astype(pd.UInt64Dtype())
    if "gene_count" in combined_qc.columns:
        combined_qc["gene_count"] = combined_qc["gene_count"].astype(pd.UInt64Dtype())
    return _reorder_qc_columns(combined_qc)
