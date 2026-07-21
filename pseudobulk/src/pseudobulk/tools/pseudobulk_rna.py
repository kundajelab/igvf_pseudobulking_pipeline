import logging
from collections import defaultdict
from collections.abc import Iterable
from concurrent.futures import (
    ThreadPoolExecutor,
    as_completed,
)
from multiprocessing import cpu_count
from pathlib import Path
from threading import Lock

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from pseudobulk import utils
from pseudobulk.types import PseudobulkName


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
    adata = ad.read_h5ad(f"{h5ad_path}")
    adata.obs["analysis_set_accession"] = h5ad_path.name.split(".", 1)[0]
    adata.obs["barcode_sample"] = adata.obs.index
    # Compute QC for each file (analysis_set_accession)
    adata.var["gene_symbol"] = adata.var.index.map(gene_ref["gene_name"])
    adata.var["mt"] = adata.var.index.map(gene_ref["mt"]).fillna(False).astype(bool)
    adata.var["ribo"] = adata.var.index.map(gene_ref["ribo"]).fillna(False).astype(bool)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def _save_obs(adata: ad.AnnData, metadata_df_x: pd.DataFrame, qc_report_path: Path) -> None:
    """Save QC info for AnnData corresponding to a raw RNA h5ad."""
    obs = pd.DataFrame(adata.obs).rename(
        columns={
            "total_counts": "rna_read_count",
            "n_genes_by_counts": "gene_count",
            "pct_counts_mt": "pct_mito",
            "pct_counts_ribo": "pct_ribo",
        },
    )  # Rename
    obs["rna_read_count"] = obs["rna_read_count"].astype(int)
    obs["annotated"] = obs["barcode_sample"].isin(set(metadata_df_x["barcode_sample"]))
    barcodes_to_pseudobulks = utils.map_barcodes_to_pseudobulks(metadata_df_x)
    obs["pseudobulk_id"] = obs["barcode_sample"].map(
        lambda barcode: barcodes_to_pseudobulks.get(barcode, "null")
    )
    # restrict to desired columns
    wanted_cols = list(utils.RNA_QC_COLUMNS)
    obs = obs[wanted_cols]
    # update the ann_data with the new obs. Unclear if this step is necessary or helpful
    adata.obs = obs
    # Save QC for all cells in analysis accession
    qc_report_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if qc_report_path.suffix == ".tsv" else ","
    obs.to_csv(f"{qc_report_path}", sep=sep, index=False)


def _save_cell_name_to_annotation_mapping(metadata_df, output_dir: Path) -> None:
    """Get mapping from cell name to annotation and save to TSV."""
    cell_name_to_annotation_df = metadata_df[
        ["pseudobulk_id", "cell_name", "annotation", "CL_id", "cell_description", "subsample"]
    ].drop_duplicates()
    cell_name_to_annotation_df.to_csv(
        f"{output_dir}/cell_name_to_annotation_mapping.tsv", sep="\t", index=False
    )


def _load_and_qc_h5ad(
    h5ad_path: Path,
    metadata_df: pd.DataFrame,
    gene_ref: pd.DataFrame,
    rna_qc_reports_dir: Path,
    logger: logging.Logger,
    log_lock: Lock,
) -> list[tuple[PseudobulkName, ad.AnnData]]:
    accession = h5ad_path.name.split(".", 1)[0]
    with log_lock:
        logger.info(f"Processing h5ad file: {h5ad_path} for accession {accession}")
    metadata_df_x: pd.DataFrame = metadata_df.loc[
        metadata_df["analysis_set_accession"] == accession, :
    ]
    adata = _load_ann_data(h5ad_path, gene_ref=gene_ref)

    qc_report_path = rna_qc_reports_dir / f"{accession}.scRNA_all_cells_QC_metrics.tsv"
    _save_obs(adata, metadata_df_x=metadata_df_x, qc_report_path=qc_report_path)

    # Get annotated data for each pseudobulk ID
    pseudobulk_adatas = []
    for pseudobulk_id, metadata_df_xcs in metadata_df_x.groupby(
        "pseudobulk_id", sort=False, group_keys=False, as_index=False
    ):
        logger.info(f"\tIn {h5ad_path}: separating out AnnData for pseudobulk_id {pseudobulk_id}")
        barcodes_xcs = set(metadata_df_xcs["barcode_sample"])
        adata_xcs = adata[adata.obs_names.isin(barcodes_xcs), :].copy()
        adata_xcs.obs["pseudobulk_id"] = pseudobulk_id
        pseudobulk_adatas.append((PseudobulkName(f"{pseudobulk_id}"), adata_xcs))
    return pseudobulk_adatas


def _load_and_qc_h5ads(
    executor: ThreadPoolExecutor,
    metadata_loc: Path,
    h5ad_paths: Iterable[Path],
    gene_ref: pd.DataFrame,
    output_dir: Path,
    rna_qc_reports_dir: Path,
    logger: logging.Logger,
    log_lock: Lock,
) -> defaultdict[PseudobulkName, list[ad.AnnData]]:
    """Load raw RNA h5ads and separate into pseudobulked AnnData objects.

    Also save observation and QC data for each accession ID into rna_qc_reports_dir.

    Args:
        metadata_loc: Path to metadata file. The file should be a tab-separated values (TSV) file
            with columns "analysis_set_accession", "barcode_sample", "subsample", "subsample_orig", and
            "annotation"
        h5ad_paths: Iterable of paths to raw RNA h5ad files.
        gene_ref: pandas DataFrame with gene info specific to the requested species
        output_dir: Path to folder to save outputs. It will have two sub-folders:
            "rna_qc_reports" will contain QC CSVs
            "pseudobulks" will contain pseudobulked h5ad and TSVs
        rna_qc_reports_dir: Path to folder to save RNA qc reports
        logger: logger object
    Returns:
        defaultdict with keys being pseudobulk IDs, and values being a list of AnnData objects with
            raw RNA data corresponding to that pseudobulk ID
    """
    # Load metadata
    metadata_df = utils.load_metadata(metadata_loc)

    # write map of cell name to annotation
    _save_cell_name_to_annotation_mapping(metadata_df=metadata_df, output_dir=output_dir)

    # Iterate through h5ads
    futures = [
        executor.submit(
            _load_and_qc_h5ad,
            h5ad_path=h5ad_path,
            metadata_df=metadata_df,
            gene_ref=gene_ref,
            rna_qc_reports_dir=rna_qc_reports_dir,
            logger=logger,
            log_lock=log_lock,
        )
        for h5ad_path in h5ad_paths
    ]
    pseudobulk_adatas: defaultdict[PseudobulkName, list[ad.AnnData]] = defaultdict(list)
    for adatas_result in as_completed(futures):
        for pseudobulk_id, adata in adatas_result.result():
            pseudobulk_adatas[pseudobulk_id].append(adata)

    return pseudobulk_adatas


def _aggregate_pseudobulk(
    pseudobulk_id: PseudobulkName,
    adatas: list[ad.AnnData],
    gene_ref: pd.DataFrame,
    rna_qc_reports_dir: Path,
    pseudobulked_rna_dir: Path,
    logger: logging.Logger,
    log_lock: Lock,
) -> None:
    """Aggregate AnnDatas for this pseudobulk ID, save pseudobulked data, counts, and QC.

    Args:
        pseudobulk_id: ID for this pseudobulk
        adatas: list of AnnData objects that correspond to this pseudobulk
        gene_ref: pandas DataFrame with gene info specific to the requested species
        rna_qc_reports_dir: Path to folder to save pseudobulked RNA QC reports
        pseudobulked_rna_dir: Path to folder to save pseudobulked RNA
        logger: logger object
    """
    with log_lock:
        logger.info(f"Aggregating pseudobulk_id: {pseudobulk_id}")
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


def pseudobulk_rna(
    *,
    input: Path,
    output_dir: Path,
    metadata_loc: Path,
    gene_info: Path,
    num_workers: int = -1,
) -> None:
    """Separate RNA h5ad files by pseudorep.

    Args:
        input: Either a folder, in which case raw input RNA files will be globbed from this folder,
            or a file, with one input RNA file per line
        output_dir: Path to folder to save outputs. It will have two sub-folders:
            "rna_qc_reports" will contain QC CSVs
            "pseudobulks" will contain pseudobulked h5ad and TSVs
        metadata_loc: Input annotations metadata file path
        at_annotation_level: Whether to split pseudobulks only at the annotation level (if True) or
            also at the subsample level (if False).
        gene_info: Path to species-specific CSV of gene info
        num_workers: Number of parallel workers to use. If <=0, use all available cores.
    """
    logger = logging.getLogger(name="pseudobulk-rna")
    num_workers = num_workers if num_workers > 0 else cpu_count()
    # Load gene information
    gene_ref = utils.read_csv(gene_info, index_col=0)
    rna_qc_reports_dir = output_dir / "rna_qc_reports"
    log_lock = Lock()
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        logger.info(f"Loading and QC-ing h5ads with {num_workers} workers.")
        # Load raw RNA h5ads, save per-accession QC info, and group by pseudobulk ID
        pseudobulk_adatas = _load_and_qc_h5ads(
            executor=executor,
            metadata_loc=metadata_loc,
            h5ad_paths=_get_h5ad_paths(input, logger=logger),
            gene_ref=gene_ref,
            output_dir=output_dir,
            rna_qc_reports_dir=rna_qc_reports_dir,
            logger=logger,
            log_lock=log_lock,
        )

        # aggregate across pseudobulks and save
        logger.info(f"Aggregating {len(pseudobulk_adatas)} pseudobulks with {num_workers} workers.")
        pseudobulked_rna_dir = output_dir / "pseudobulks"
        futures = [
            executor.submit(
                _aggregate_pseudobulk,
                pseudobulk_id=pseudobulk_id,
                adatas=adatas,
                gene_ref=gene_ref,
                rna_qc_reports_dir=rna_qc_reports_dir,
                pseudobulked_rna_dir=pseudobulked_rna_dir,
                logger=logger,
                log_lock=log_lock,
            )
            for pseudobulk_id, adatas in pseudobulk_adatas.items()
        ]
        for future in as_completed(futures):
            future.result()  # raise any exceptions in worker threads
