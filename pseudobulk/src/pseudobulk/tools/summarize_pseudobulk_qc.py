import logging
from pathlib import Path

import numpy as np
import pandas as pd

from pseudobulk import utils
from pseudobulk.types import (
    COUNTS_MATRIX,
    FailureAction,
    FailureHandler,
    PseudobulkName,
)


def _load_and_subset_metadata(metadata_loc: Path, pseudobulk: PseudobulkName) -> pd.DataFrame:
    """Load metadata and subset to the specified pseudobulk."""
    metadata_df = utils.load_metadata(metadata_loc)
    return metadata_df.loc[metadata_df["pseudobulk_id"] == pseudobulk, :]


def _load_and_subset_atac_qc(
    atac_qc_dir: Path,
    metadata_pseudobulk_df: pd.DataFrame,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, COUNTS_MATRIX]:
    """Load ATAC QC and subset to the specified pseudobulk."""
    pseudobulk_id = metadata_pseudobulk_df["pseudobulk_id"].iloc[0]

    def _row_filter(df: pd.DataFrame) -> pd.Series[bool]:
        return df["pseudobulk_id"] == pseudobulk_id

    pseudobulk_atac_qc, pseudobulk_tss_matrix = utils.load_atac_qc(
        atac_qc_dir=atac_qc_dir, need_tss=True, logger=logger, row_filter=_row_filter
    )
    subsample = metadata_pseudobulk_df["subsample"].iloc[0]
    pseudobulk_atac_qc["subsample"] = subsample
    return pseudobulk_atac_qc, pseudobulk_tss_matrix


def _compute_pseudobulk_combined_qc(
    pseudobulk: PseudobulkName,
    pseudobulk_atac_qc: pd.DataFrame,
    rna_qc: Path | None,
    frip_per_cell: Path | None,
    failure_handler: FailureHandler,
    logger: logging.Logger,
) -> pd.DataFrame:
    logger.info(f"Computing pseudobulk_combined_qc for {pseudobulk}")
    logger.info(f"Loading RNA QC for {pseudobulk}")
    pseudobulk_rna_qc = (
        utils.read_csv(rna_qc)
        if rna_qc is not None and rna_qc.exists()
        else pd.DataFrame([], columns=pd.Index(utils.RNA_QC_COLUMNS))
    )
    # copy pseudobulk_atac_qc so that changes aren't propagated to the caller
    pseudobulk_atac_qc = pseudobulk_atac_qc.copy(deep=True)
    # NOTE: remove raw columns thatare only used for pseudobulk-lvl QC
    pseudobulk_atac_qc = pseudobulk_atac_qc.loc[
        :, [x for x in pseudobulk_atac_qc.columns if not x.startswith("raw-")]
    ]
    if len(pseudobulk_atac_qc) > 0 and frip_per_cell is not None and frip_per_cell.exists():
        pseudobulk_frip_qc = utils.read_csv(frip_per_cell, names=["barcode_sample", "frip"])
        # Combine ATAC QC with frip data
        pseudobulk_atac_qc = pd.merge(
            pseudobulk_atac_qc, pseudobulk_frip_qc, how="outer", on="barcode_sample"
        )
    pseudobulk_combined_qc = utils.merge_rna_and_atac_qc(
        identifier=pseudobulk,
        rna_qc=pseudobulk_rna_qc,
        atac_qc=pseudobulk_atac_qc,
        logger=logger,
    )

    # Confirm that ATAC and RNA cells match
    if (
        len(pseudobulk_rna_qc) > 0
        and len(pseudobulk_atac_qc) > 0
        and not (len(pseudobulk_rna_qc) == len(pseudobulk_atac_qc) == len(pseudobulk_combined_qc))
    ):
        failure_handler.handle_failure(
            f"ATAC and RNA pseudobulk {pseudobulk} cell sets do not match!"
        )

    # return subset of columns
    combined_cols = [
        "analysis_set_accession",
        "barcode_sample",
        "subsample",
        "rna_read_count",
        "gene_count",
        "pct_mito",
        "pct_ribo",
        "num_frags",
        "pct_duplicated_reads",
        "nucleosomal_signal",
        "tss_enrichment",
        "frip",
    ]
    return pseudobulk_combined_qc.loc[:, combined_cols]


def _compute_pseudobulk_qc_summary(
    pseudobulk: PseudobulkName,
    pseudobulk_atac_qc: pd.DataFrame,
    pseudobulk_tss_matrix: COUNTS_MATRIX,
    cell_name: str,
    subsample: str | None,
    pseudobulk_combined_qc: pd.DataFrame,
    pseudobulk_counts: Path | None,
    fragments_per_cell: Path | None,
    fragments_in_peaks_per_cell: Path | None,
    logger: logging.Logger,
) -> pd.Series:
    # Compute pseudobulk QC summary
    logger.info(f"Computing pseudobulk_qc_summary for {pseudobulk}")
    pseudobulk_qc_summary = pd.Series()
    pseudobulk_qc_summary["directory_name"] = pseudobulk
    pseudobulk_qc_summary["pseudobulk"] = f"{cell_name}-{subsample}"
    pseudobulk_qc_summary["cell_name"] = cell_name
    pseudobulk_qc_summary["subsample"] = subsample
    pseudobulk_qc_summary["num_cells"] = pseudobulk_combined_qc.shape[0]
    # RNA
    if pseudobulk_counts is None or not pseudobulk_counts.exists():
        pseudobulk_qc_summary["rna_read_count"] = None
        pseudobulk_qc_summary["gene_count"] = None
        pseudobulk_qc_summary["pct_mito"] = float("nan")
        pseudobulk_qc_summary["pct_ribo"] = float("nan")
    else:
        pseudobulk_rna_exp = utils.read_csv(pseudobulk_counts)
        pseudobulk_qc_summary["rna_read_count"] = pseudobulk_rna_exp["counts"].sum()
        pseudobulk_qc_summary["gene_count"] = np.count_nonzero(pseudobulk_rna_exp["counts"])
        pseudobulk_qc_summary["pct_mito"] = (
            pseudobulk_rna_exp[pseudobulk_rna_exp["mt"]]["counts"].sum()
            / pseudobulk_qc_summary["rna_read_count"]
        ) * 100
        pseudobulk_qc_summary["pct_ribo"] = (
            pseudobulk_rna_exp[pseudobulk_rna_exp["ribo"]]["counts"].sum()
            / pseudobulk_qc_summary["rna_read_count"]
        ) * 100
    # ATAC
    pseudobulk_qc_summary["num_frags"] = pseudobulk_atac_qc["num_frags"].sum()
    pseudobulk_qc_summary["pct_duplicated_reads"] = (
        pseudobulk_atac_qc["raw-num_dup_reads"].sum() / pseudobulk_atac_qc["raw-num_reads"].sum()
    ) * 100
    pseudobulk_qc_summary["nucleosomal_signal"] = (
        1 + pseudobulk_atac_qc["raw-mono_nucleosomal_frags"].sum()
    ) / (1 + pseudobulk_atac_qc["raw-nucleosome_free_frags"].sum())
    # ATAC - TSS enrichment
    TSS_half_window = 2000
    TSS_half_smooth_window = 5
    pseudobulk_tss_insertions = pseudobulk_tss_matrix.sum(axis=0, dtype=np.uint64)
    tss_insertions_flank_mean = (
        pseudobulk_tss_insertions[:100].sum() + pseudobulk_tss_insertions[-100:].sum()
    ) / 200
    tss_insertions_center = pseudobulk_tss_insertions[
        TSS_half_window - TSS_half_smooth_window : TSS_half_window + TSS_half_smooth_window + 1
    ].mean()
    pseudobulk_qc_summary["tss_enrichment"] = tss_insertions_center / (
        tss_insertions_flank_mean + 0.1
    )  # add 0.1 like snapatac2 to avoid division by zero
    # ATAC - FRIP
    if (
        fragments_per_cell is None
        or not fragments_per_cell.exists()
        or fragments_in_peaks_per_cell is None
        or not fragments_in_peaks_per_cell.exists()
    ):
        pseudobulk_qc_summary["frip"] = None
    else:
        fragments_per_cell_df = utils.read_csv(
            fragments_per_cell, names=["barcode_sample", "num_fragments"]
        )
        fragments_in_peaks_per_cell_df = utils.read_csv(
            fragments_in_peaks_per_cell,
            names=["barcode_sample", "num_fragments_in_peaks"],
        )
        merged_frip = pd.merge(
            fragments_per_cell_df,
            fragments_in_peaks_per_cell_df,
            how="left",
            on="barcode_sample",
        )
        merged_frip["num_fragments_in_peaks"] = merged_frip["num_fragments_in_peaks"].fillna(0)
        pseudobulk_qc_summary["frip"] = (
            merged_frip["num_fragments_in_peaks"].sum() / merged_frip["num_fragments"].sum()
        )
    return pseudobulk_qc_summary


def summarize_pseudobulk_qc(
    *,
    pseudobulk: str,
    metadata_loc: Path,
    atac_qc_dir: Path,
    pseudobulk_qc_out: Path,
    qc_summary_out: Path,
    rna_qc: Path | None = None,
    pseudobulk_counts: Path | None = None,
    frip_per_cell: Path | None = None,
    fragments_per_cell: Path | None = None,
    fragments_in_peaks_per_cell: Path | None = None,
    failure_actions: list[FailureAction] = [FailureAction.warning, FailureAction.sentinal],
) -> None:
    """Summarize ATAC and RNA QC for the specified pseudobulk.

    Produces
    -pseudobulk_combined_qc, which is the final per-pseudobulk QC file
    -pseudobulk_qc_summary, which will (in a separate step) be concatenated into an overall QC file

    Args:
        pseudobulk: ID of the pseudobulk to process
        metadata_loc: Path to metadata file. The file should be a tab-separated values (TSV) file
            with columns "analysis_set_accession", "barcode_sample", "subsample", "subsample_orig", and
            "annotation"
        atac_qc_dir: Path to folder with ATAC QC TSVs and .npy files with Transcription Start Sites,
            produced by the split-fragments tool.
        rna_qc: Path to RNA QC TSV produced by pseudobulk-rna tool.
        pseudobulk_counts: Path to RNA counts TSV produced by pseudobulk-rna tool.
        frip_per_cell: Path to TXT file produced by CALL_PEAKS module.
        fragments_per_cell: Path to TXT file produced by CALL_PEAKS module.
        fragments_in_peaks_per_cell: Path to TXT file produced by CALL_PEAKS module.
        output_dir: Path,
        failure_actions: List of potential actions to take if data integrity checks fail.
    """
    logger = logging.getLogger(name=f"{__package__} aggregate-pseudobulk-qc")
    failure_handler = FailureHandler(
        actions=failure_actions, sentinal_file_name=Path("ATAC_RNA_MISMATCH.log"), logger=logger
    )
    pseudobulk = PseudobulkName(pseudobulk)
    logger.info(f"Loading metadata for {pseudobulk}")
    metadata_pseudobulk_df = _load_and_subset_metadata(
        metadata_loc=metadata_loc, pseudobulk=pseudobulk
    )
    logger.info(f"Loading ATAC QC for {pseudobulk}")
    pseudobulk_atac_qc, pseudobulk_tss_matrix = _load_and_subset_atac_qc(
        atac_qc_dir=atac_qc_dir,
        metadata_pseudobulk_df=metadata_pseudobulk_df,
        logger=logger,
    )
    pseudobulk_combined_qc = _compute_pseudobulk_combined_qc(
        pseudobulk=pseudobulk,
        pseudobulk_atac_qc=pseudobulk_atac_qc,
        rna_qc=rna_qc,
        frip_per_cell=frip_per_cell,
        failure_handler=failure_handler,
        logger=logger,
    )
    pseudobulk_combined_qc.to_csv(f"{pseudobulk_qc_out}", sep="\t", index=False)

    cell_name = metadata_pseudobulk_df["cell_name"].iloc[0]
    subsample = metadata_pseudobulk_df["subsample"].iloc[0]
    pseudobulk_qc_summary = _compute_pseudobulk_qc_summary(
        pseudobulk=pseudobulk,
        pseudobulk_atac_qc=pseudobulk_atac_qc,
        pseudobulk_tss_matrix=pseudobulk_tss_matrix,
        cell_name=cell_name,
        subsample=subsample,
        pseudobulk_combined_qc=pseudobulk_combined_qc,
        pseudobulk_counts=pseudobulk_counts,
        fragments_per_cell=fragments_per_cell,
        fragments_in_peaks_per_cell=fragments_in_peaks_per_cell,
        logger=logger,
    )
    # write as row of TSV, dropping index
    (
        pd.DataFrame(pseudobulk_qc_summary)
        .transpose()
        .to_csv(f"{qc_summary_out}", sep="\t", index=False)
    )
