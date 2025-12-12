import argparse
import os
import shutil
import subprocess

import numpy as np
import pandas as pd

from utils import load_metadata


##################
# MAIN FUNCTIONS #
##################
def aggregate_files(data_dir, metadata_loc, at_annotation_level):
    # Load metadata
    metadata_df = load_metadata(metadata_loc)
    print(metadata_df)
    # Load combined atac qc (generated per analysis accession, not by pseudobulk)
    atac_combined_qc_list = []
    tss_list = []
    for x in [x for x in os.listdir(f"{data_dir}/atac_qc_reports") if x.endswith(".tsv")]:
        analysis_accession = x.split(".")[0]
        atac_combined_qc_list.append(pd.read_csv(f"{data_dir}/atac_qc_reports/{x}", sep="\t"))
        tss_list.append(np.load(f"{data_dir}/atac_qc_reports/{analysis_accession}_tss_matrix.npy"))
    atac_combined_qc = pd.concat(atac_combined_qc_list, axis=0)
    print(atac_combined_qc)
    tss_matrix = np.vstack(tss_list)
    # Find all pseudobulks (use pseudobulked RNA files as proxy to find pseudobulks)
    pseudobulks = [x.split(".")[0] for x in os.listdir(f"{data_dir}/pseudobulked_rna") if x.endswith(".h5ad")]
    # Process each pseudobulk
    pseudobulk_qc_rows = []
    for pseudobulk in pseudobulks:
        print(f"Processing pseudobulk: {pseudobulk}")
        if at_annotation_level:
            annotation_type, annotation_value = tuple(pseudobulk.split("-"))
            metadata_pseudobulk_df = metadata_df[metadata_df[annotation_type] == annotation_value]
        else:
            annotation_type, annotation_value, subsample = tuple(pseudobulk.split("-"))
            metadata_pseudobulk_df = metadata_df[(metadata_df[annotation_type] == annotation_value) & (metadata_df["subsample"] == subsample)]
        # Create pseudobulk directory
        os.makedirs(f"{data_dir}/pseudobulks/{pseudobulk}", exist_ok=True)
        # Move ATAC files
        subprocess.run([f"gzip -c {data_dir}/pseudobulked_fragments/{pseudobulk}-sorted.tsv > {data_dir}/pseudobulks/{pseudobulk}/fragments.tsv.gz"], shell=True)
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-raw_insertions.bw", f"{data_dir}/pseudobulks/{pseudobulk}/raw_insertions.bw")
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-peaks_overlap_filtered.narrowPeak", f"{data_dir}/pseudobulks/{pseudobulk}/peaks.narrowPeak")
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-pval.bw", f"{data_dir}/pseudobulks/{pseudobulk}/peaks_minuslog10pval.bw")
        # Move RNA files
        shutil.copy(f"{data_dir}/pseudobulked_rna/{pseudobulk}.h5ad", f"{data_dir}/pseudobulks/{pseudobulk}/rna_counts_mtx.h5ad")
        shutil.copy(f"{data_dir}/pseudobulked_rna/{pseudobulk}-pseudobulked_counts.tsv", f"{data_dir}/pseudobulks/{pseudobulk}/pseudobulk_expression.tsv")
        # RNA QC - load
        pseudobulk_rna_qc = pd.read_csv(f"{data_dir}/rna_qc_reports/{pseudobulk}-pseudobulked_cell_QC_metrics.csv", sep="\t")
        # ATAC QC - build from combined atac qc
        pseudobulk_atac_qc_chunks = []
        pseudobulk_tss_chunks = []
        for accession in metadata_pseudobulk_df['analysis_accession'].unique().tolist():
            metadata_pseudobulk_accession = metadata_pseudobulk_df[metadata_pseudobulk_df['analysis_accession'] == accession]
            atac_combined_qc_chunk = atac_combined_qc[(atac_combined_qc['analysis_accession'] == accession) & (atac_combined_qc['barcode'].isin(set(metadata_pseudobulk_accession['barcode'])))].copy()
            barcode_to_subsample = dict(zip(metadata_pseudobulk_accession['barcode'], metadata_pseudobulk_accession['subsample']))
            atac_combined_qc_chunk["subsample"] = [barcode_to_subsample[x] for x in atac_combined_qc_chunk['barcode']]
            pseudobulk_atac_qc_chunks.append(atac_combined_qc_chunk)
            pseudobulk_tss_chunks.append(tss_matrix[atac_combined_qc_chunk.index])
        pseudobulk_atac_qc = pd.concat(pseudobulk_atac_qc_chunks, axis=0) # NOTE: contains 'raw-' columns
        pseudobulk_tss_matrix = np.vstack(pseudobulk_tss_chunks)
        pseudobulk_atac_qc_raw = pseudobulk_atac_qc.copy()
        pseudobulk_atac_qc = pseudobulk_atac_qc[[x for x in pseudobulk_atac_qc.columns if not x.startswith("raw-")]] # NOTE: raw columns just used for pseudobulk-lvl QC
        pseudobulk_frip_qc = pd.read_csv(f"{data_dir}/peaks/{pseudobulk}-frip_per_cell.txt", sep=" ", names=["barcode", "frip"])
        # Combine QC
        pseudobulk_atac_qc = pd.merge(pseudobulk_atac_qc, pseudobulk_frip_qc, how="outer", on="barcode")
        pseudobulk_combined_qc = pd.merge(pseudobulk_rna_qc, pseudobulk_atac_qc, how="outer", on=["analysis_accession", "barcode", "annotated"])
        pseudobulk_combined_qc = pseudobulk_combined_qc[["analysis_accession", "barcode", "subsample", "rna_read_count", "gene_count", "pct_mito", "pct_ribo", "num_frags", "pct_duplicated_reads", "nucleosomal_signal", "tss_enrichment", "frip"]].copy()
        print(pseudobulk_combined_qc)
        pseudobulk_combined_qc.to_csv(f"{data_dir}/pseudobulks/{pseudobulk}/per_cell_qc.tsv", sep="\t", index=False)
        # Confirm that ATAC and RNA cells match
        if not (len(pseudobulk_rna_qc) == len(pseudobulk_atac_qc) == len(pseudobulk_combined_qc)):
            print("!!!!! WARNING: ATAC and RNA pseudobulk cell sets do not match !!!!!")
            with open(f"{data_dir}/pseudobulks/{pseudobulk}/ATAC_RNA_MISMATCH.log", 'w') as f:
                f.write(f"WARNING: ATAC {len(pseudobulk_combined_qc[pseudobulk_combined_qc['found_in_atac']])} and RNA {len(pseudobulk_combined_qc[pseudobulk_combined_qc['found_in_rna']])} pseudobulk cell sets do not match!\n")
        # Compute pseudobulk QC summary
        pseudobulk_qc_summary = dict()
        pseudobulk_qc_summary["pseudobulk"] = pseudobulk
        pseudobulk_qc_summary["annotation_level"] = annotation_type
        pseudobulk_qc_summary["annotation"] = annotation_value
        if not at_annotation_level:
            pseudobulk_qc_summary["subsample"] = subsample
        pseudobulk_qc_summary["num_cells"] = pseudobulk_combined_qc.shape[0]
        # RNA
        pseudobulk_rna_exp = pd.read_csv(f"{data_dir}/pseudobulks/{pseudobulk}/pseudobulk_expression.tsv", sep="\t")
        pseudobulk_qc_summary["rna_read_count"] = pseudobulk_rna_exp["counts"].sum()
        pseudobulk_qc_summary["gene_count"] = np.count_nonzero(pseudobulk_rna_exp["counts"])
        pseudobulk_qc_summary["pct_mito"] = (pseudobulk_rna_exp[pseudobulk_rna_exp["mt"]]["counts"].sum() / pseudobulk_qc_summary["rna_read_count"]) * 100
        pseudobulk_qc_summary["pct_ribo"] = (pseudobulk_rna_exp[pseudobulk_rna_exp["ribo"]]["counts"].sum() / pseudobulk_qc_summary["rna_read_count"]) * 100
        # ATAC
        pseudobulk_qc_summary["num_frags"] = pseudobulk_atac_qc_raw["num_frags"].sum()
        pseudobulk_qc_summary["pct_duplicated_reads"] = (pseudobulk_atac_qc_raw["raw-num_dup_reads"].sum() / pseudobulk_atac_qc_raw["raw-num_reads"].sum()) * 100
        pseudobulk_qc_summary["nucleosomal_signal"] = (1 + pseudobulk_atac_qc_raw["raw-mono_nucleosomal_frags"].sum()) / (1 + pseudobulk_atac_qc_raw["raw-nucleosome_free_frags"].sum())
        pseudobulk_tss_insertions = pseudobulk_tss_matrix.sum(axis=0)
        tss_insertions_flank_mean = (np.sum(pseudobulk_tss_insertions[:100]) + np.sum(pseudobulk_tss_insertions[-100:])) / 200
        tss_insertions_center = np.mean(pseudobulk_tss_insertions[1000-2:1000+3])
        pseudobulk_qc_summary["tss_enrichment"] = (tss_insertions_center / (tss_insertions_flank_mean if tss_insertions_flank_mean > 0 else 1))
        # ATAC - FRIP
        fragments_per_cell_df = pd.read_csv(f"{data_dir}/peaks/{pseudobulk}-fragments_per_cell.txt", sep=" ", names=["barcode", "num_fragments"])
        fragments_in_peaks_per_cell_df = pd.read_csv(f"{data_dir}/peaks/{pseudobulk}-fragments_in_peaks_per_cell.txt", sep=" ", names=["barcode", "num_fragments_in_peaks"])
        assert(len(fragments_per_cell_df) == len(fragments_in_peaks_per_cell_df))
        merged_frip = pd.merge(fragments_per_cell_df, fragments_in_peaks_per_cell_df, how="inner", on="barcode")
        pseudobulk_qc_summary["frip"] = np.sum(merged_frip["num_fragments_in_peaks"]) / np.sum(merged_frip["num_fragments"])
        # Append to list
        pseudobulk_qc_rows.append(pseudobulk_qc_summary)
    # Save pseudobulk QC
    pseudobulk_qc_df = pd.DataFrame(pseudobulk_qc_rows)
    print(pseudobulk_qc_df)
    pseudobulk_qc_df.to_csv(f"{data_dir}/pseudobulk_qc.tsv", sep="\t", index=False)
    # Combine QC per analysis accession
    atac_combined_qc = atac_combined_qc[[x for x in atac_combined_qc.columns if not x.startswith("raw-")]] # NOTE: don't need raw columns after pseudobulk processing
    for accession in metadata_df['analysis_accession'].unique().tolist():
        print(f"Processing analysis accession: {accession}")
        atac_qc = atac_combined_qc[atac_combined_qc['analysis_accession'] == accession].copy()
        atac_qc["found_in_atac"] = True
        print(atac_qc)
        rna_qc = pd.read_csv(f"{data_dir}/rna_qc_reports/{accession}-scRNA_all_cells_QC_metrics.csv", sep="\t")
        rna_qc["found_in_rna"] = True
        print(rna_qc)
        # Combine QC
        combined_qc = pd.merge(rna_qc, atac_qc, how="outer", on=["analysis_accession", "barcode", "annotated"])
        combined_qc["found_in_rna"] = combined_qc["found_in_rna"].fillna(False).astype(bool)
        combined_qc["found_in_atac"] = combined_qc["found_in_atac"].fillna(False).astype(bool)
        combined_qc = combined_qc[["analysis_accession", "barcode", "annotated", "found_in_rna", "found_in_atac"] + [x for x in combined_qc.columns if x not in ["analysis_accession", "barcode", "annotated", "found_in_rna", "found_in_atac"]]]
        combined_qc.to_csv(f"{data_dir}/analysis_accession_qc_reports/{accession}_per_cell_qc.tsv", sep="\t", index=False)
    # CONFIRM COMPLETION
    with open(f"{data_dir}/step5_complete.txt", 'w') as f:
        f.write("Step 5: Aggregate Outputs completed successfully!\n")


def main():
    parser = argparse.ArgumentParser(description="aggregate output files")

    # Define the expected flags
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')
    parser.add_argument('-a', '--at_annotation_level', type=str, required=True, help='At annotation level flag')

    # Parse the arguments
    args = parser.parse_args()

    if args.at_annotation_level == "True":
        args.at_annotation_level = True
    elif args.at_annotation_level == "False":
        args.at_annotation_level = False
    else:
        raise ValueError("at_annotation_level must be 'True' or 'False'")

    # Process file
    aggregate_files(args.datadir, args.metadata, args.at_annotation_level)
    

if __name__ == "__main__":
    main()

