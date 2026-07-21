from collections.abc import (
    Hashable,
    Mapping,
    Sequence,
)
from contextlib import contextmanager
from contextlib import nullcontext
from logging import Logger
from threading import Lock
from pathlib import Path
from types import MappingProxyType
from typing import (
    Callable,
    Final,
    Generator,
    Literal,
    TextIO,
    overload,
)

import numpy as np
import pandas as pd
import scipy.sparse

from pseudobulk.barcode_qc import BarcodeQc
from pseudobulk.types import (
    COUNTS_ARRAY,
    COUNTS_MATRIX,
    POS_ARRAY,
    POS_DTYPE,
    Barcode,
    Contig,
    PseudobulkName,
)
from pseudobulk.tss import update_insertions_range


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

FRONT_COLUMNS: Final[tuple[str, ...]] = (
    "analysis_set_accession",
    "barcode_sample",
    "annotated",
    "found_in_rna",
    "found_in_atac",
    "pseudobulk_id",
)
"""Columns that should come first in combined RNA + ATAC QC."""


def load_metadata(metadata_loc: Path) -> pd.DataFrame:
    """Load metadata file and perform basic checks."""
    mandatory_columns = {
        "barcode_sample",
        "cell_name",
        "cell_description",
        "CL_id",
        "CL_term_name",
        "subsample",
        "analysis_set_accession",
    }
    df = read_csv(metadata_loc)
    for col in mandatory_columns:
        if col not in df.columns:
            raise ValueError(f"metadata must contain '{col}' column")
    for subsample in df["subsample"].unique():
        if "-" in subsample:
            raise ValueError(
                f"'subsample' column contained '{subsample}' containing invalid annotation values"
            )
    _map_annotations_and_pseudobulk_ids(df)
    return df


def _map_annotations_and_pseudobulk_ids(metadata_df) -> None:
    """
    Update metadata_df in place with new 'annotation' and 'pseudobulk_id' columns.
    """
    # Cell name to annotation mapping
    cell_name_to_annotation_dict = {
        x: f"annotation_{i}" for i, x in enumerate(sorted(metadata_df["cell_name"].unique()))
    }
    # Map cell names to annotations in metadata_df
    metadata_df["annotation"] = metadata_df["cell_name"].map(cell_name_to_annotation_dict)
    metadata_df["pseudobulk_id"] = metadata_df.apply(
        lambda row: f"{row['annotation']}-{row['subsample']}", axis=1
    )


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
        # for idx in range(tss_left_idx, tss_right_idx):
        #     tss = tss_vec[idx]
        #     strand = strand_vec[idx]
        #     tss_insertions[half_window + (position - tss) * strand] += 1
        update_insertions_range(
            tss_insertions, tss_vec, strand_vec, tss_left_idx, tss_right_idx, half_window, position
        )

        # this would be the vectorized version, but it is actually slower, likely because of one of
        # these issues:
        # 1. numpy tries to vectorize over threads and creates contention / overhead
        # 2. intermediate arrays are created and numpy malloc blocks across threads.
        # When python 3.15 is stable, it may be worth revisiting, although it may require more than
        # one "index_buffer" to prevent creation of intemediate arrays
        # index_buffer[:tss_idx_delta] = offset + strand_vec[tss_left_idx:tss_right_idx] * (
        #     position - tss_vec[tss_left_idx:tss_right_idx]
        # )
        # tss_insertions[index_buffer[:tss_idx_delta]] += 1
        # Rather than vectorizing, it may also be possible to create a ufunc to do this in cython


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


def read_csv(path: Path, sep: str | None = None, **kwargs) -> pd.DataFrame:
    if sep is None:
        match path.suffixes:
            case [*_parts, ".tsv"] | [*_parts, ".tsv", ".gz"]:
                sep = "\t"
            case [*_parts, ".csv"] | [*_parts, ".csv", ".gz"]:
                sep = ","
            case _:
                raise ValueError(
                    f"Could not infer separator for file {path} with suffix {path.suffix}. Please "
                    "specify separator with 'sep' argument."
                )
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
    log_lock: Lock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
) -> pd.DataFrame: ...


@overload
def load_atac_qc(
    atac_qc_dir: Path,
    need_tss: Literal[True],
    logger: Logger,
    log_lock: Lock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
) -> tuple[pd.DataFrame, COUNTS_MATRIX]: ...


def load_atac_qc(
    atac_qc_dir: Path,
    need_tss: bool,
    logger: Logger,
    log_lock: Lock | None = None,
    row_filter: Callable[[pd.DataFrame], pd.Series[bool]] | None = None,
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None = None,
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

        atac_combined_qc = pd.DataFrame([], pd.Index(atac_cols))
        tss_counts = scipy.sparse.csr_array(np.empty((0, 4001), np.uint16))
    else:
        atac_combined_qc = pd.concat(atac_combined_qc_list, axis=0, ignore_index=True)
        tss_counts = scipy.sparse.vstack(tss_list) if need_tss else None

    return (atac_combined_qc, tss_counts) if need_tss else atac_combined_qc


def _reorder_qc_columns(combined_qc: pd.DataFrame) -> pd.DataFrame:
    """Reorder columns to move selected columns to the front."""
    # Move shared columns to the front, and fill out the remaining columns in order afterwards
    col_order = [x for x in FRONT_COLUMNS if x in combined_qc.columns] + [
        x for x in combined_qc.columns if x not in set(FRONT_COLUMNS)
    ]
    # do the reorder
    return combined_qc[col_order]


def merge_rna_and_atac_qc(
    identifier: str,
    rna_qc: pd.DataFrame,
    atac_qc: pd.DataFrame,
    logger: Logger,
    log_lock: Lock | None = None,
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
    return _reorder_qc_columns(combined_qc)
