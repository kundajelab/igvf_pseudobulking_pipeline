import argparse
import os
import shutil

import pandas as pd

from utils import load_metadata


##################
# MAIN FUNCTIONS #
##################
def aggregate_files(data_dir, metadata_loc):
    # Load metadata
    metadata_df = load_metadata(metadata_loc)
    print(metadata_df)
    # Load combined atac qc
    atac_combined_qc = pd.concat([pd.read_csv(f"{data_dir}/atac_qc_reports/{x}", sep="\t") for x in os.listdir(f"{data_dir}/atac_qc_reports")], axis=0)
    print(atac_combined_qc)
    # Find all pseudobulks (use pseudobulked RNA files as proxy to find pseudobulks)
    pseudobulks = [x.split(".")[0] for x in os.listdir(f"{data_dir}/pseudobulked_rna") if x.endswith(".h5ad")]
    # Process each pseudobulk
    for pseudobulk in pseudobulks:
        print(f"Processing pseudobulk: {pseudobulk}")
        annotation_type, annotation_value, subsample = tuple(pseudobulk.split("-"))
        metadata_pseudobulk_df = metadata_df[(metadata_df[annotation_type] == annotation_value) & (metadata_df["subsample"] == subsample)]
        # Create pseudobulk directory
        os.mkdir(f"{data_dir}/pseudobulks/{pseudobulk}")
        # Move ATAC files
        shutil.copy(f"{data_dir}/pseudobulked_fragments/{pseudobulk}-sorted.tsv", f"{data_dir}/pseudobulks/{pseudobulk}/fragments.tsv")
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-raw_insertions.bw", f"{data_dir}/pseudobulks/{pseudobulk}/raw_insertions.bw")
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-peaks_overlap_filtered.narrowPeak", f"{data_dir}/pseudobulks/{pseudobulk}/peaks.narrowPeak")
        shutil.copy(f"{data_dir}/peaks/{pseudobulk}-pval.bw", f"{data_dir}/pseudobulks/{pseudobulk}/peaks_minuslog10pval.bw")
        # Move RNA files
        shutil.copy(f"{data_dir}/pseudobulked_rna/{pseudobulk}.h5ad", f"{data_dir}/pseudobulks/{pseudobulk}/rna_counts_mtx.h5ad")
        shutil.copy(f"{data_dir}/pseudobulked_rna/{pseudobulk}-pseudobulked_counts.tsv", f"{data_dir}/pseudobulks/{pseudobulk}/pseudobulk_expression.tsv")
        # Prepare QC file
        rna_qc = pd.read_csv(f"{data_dir}/rna_qc_reports/{pseudobulk}-pseudobulked_cell_QC_metrics.csv", sep="\t")
        print(rna_qc)
        atac_qc = []
        for accession in metadata_pseudobulk_df['analysis_accession'].unique().tolist():
            metadata_subset = metadata_pseudobulk_df[metadata_pseudobulk_df['analysis_accession'] == accession]
            atac_qc.append(atac_combined_qc[(atac_combined_qc['analysis_accession'] == accession) & (atac_combined_qc['barcode'].isin(set(metadata_subset['barcode'])))])
        atac_qc = pd.concat(atac_qc, axis=0)
        print(atac_qc)
        frip_qc = pd.read_csv(f"{data_dir}/peaks/{pseudobulk}-frip_per_cell.txt", sep=" ", names=["barcode", "frip"])
        print(frip_qc)
        # Combine QC
        atac_qc = pd.merge(atac_qc, frip_qc, how="outer", on="barcode")
        print(atac_qc)
        combined_qc = pd.merge(rna_qc, atac_qc, how="outer", on=["analysis_accession", "barcode", "annotated"])
        print(combined_qc)
        combined_qc.to_csv(f"{data_dir}/pseudobulks/{pseudobulk}/per_cell_qc.tsv", sep="\t", index=False)
    # Combine QC per analysis accession
    for accession in metadata_df['analysis_accession'].unique().tolist():
        print(f"Processing analysis accession: {accession}")
        atac_qc = atac_combined_qc[atac_combined_qc['analysis_accession'] == accession]
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
        print(combined_qc)
        combined_qc.to_csv(f"{data_dir}/analysis_accession_qc_reports/{accession}_per_cell_qc.tsv", sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="aggregate output files")

    # Define the expected flags
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    aggregate_files(args.datadir, args.metadata)
    

if __name__ == "__main__":
    main()

