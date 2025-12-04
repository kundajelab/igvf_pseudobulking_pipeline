import argparse
from collections import defaultdict
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from utils import load_metadata


##################
# MAIN FUNCTIONS #
##################
def process_h5ad(data_dir, metadata_loc, geneinfo_loc):
    # Load metadata
    metadata_df = load_metadata(metadata_loc)
    print(metadata_df)
    # Load gene information
    gene_ref = pd.read_csv(geneinfo_loc, index_col=0)
    print(gene_ref)
    # Iterate through h5ads
    pseudobulk_adatas = defaultdict(list)
    pseudobulk_qcs = defaultdict(list)
    for x in os.listdir(f"{data_dir}/raw_rna"):
        # Load h5ads
        x_name = x.split(".")[0]
        metadata_df_x = metadata_df[metadata_df["analysis_accession"] == x_name]
        adata = ad.read_h5ad(f"{data_dir}/raw_rna/{x}")
        adata.obs['analysis_accession'] = x_name
        adata.obs['barcode'] = adata.obs.index
        print(x_name)
        print(adata)
        # Compute QC for each file (analysis_accession)
        adata.var['gene_symbol'] = adata.var.index.map(gene_ref['gene_name'])
        adata.var['mt'] = adata.var.index.map(gene_ref['mt']).fillna(False).astype(bool)
        adata.var['ribo'] = adata.var.index.map(gene_ref['ribo']).fillna(False).astype(bool)
        sc.pp.calculate_qc_metrics(
            adata, 
            qc_vars=['mt', 'ribo'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        adata.obs.rename(columns={
            'total_counts': 'read_count',
            'n_genes_by_counts': 'gene_count',
            'pct_counts_mt': 'pct_mito',
            'pct_counts_ribo': 'pct_ribo'
        }, inplace=True) # Rename
        adata.obs['read_count'] = adata.obs['read_count'].astype(int)
        adata.obs['annotated'] = adata.obs['barcode'].isin(set(metadata_df_x['barcode']))
        adata.obs = adata.obs[['analysis_accession', 'barcode', 'annotated', 'read_count', 'gene_count', 'pct_mito', 'pct_ribo']]
        # Save QC for all cells in analysis accession
        adata.obs.to_csv(f"{data_dir}/rna_qc_reports/{x_name}-scRNA_all_cells_QC_metrics.csv", sep="\t", index=False)
        # Go through annotation columns --> make pseudobulks
        for c in metadata_df_x.columns:
            if not c.startswith("annotation"):
                continue
            for c_value in metadata_df_x[c].unique().tolist():
                metadata_df_xc = metadata_df_x[metadata_df_x[c] == c_value]
                for s in metadata_df_xc["subsample"].unique().tolist():
                    metadata_df_xcs = metadata_df_xc[metadata_df_xc["subsample"] == s]
                    print(f"--- {x} {c} {c_value} {s}")
                    print(metadata_df_xcs)
                    pseudobulk_name = f"{c}-{c_value}-{s}"
                    barcodes_xcs = set(metadata_df_xcs["barcode"])
                    adata_xcs = adata[adata.obs_names.isin(barcodes_xcs), :].copy()
                    pseudobulk_adatas[pseudobulk_name].append(adata_xcs)
                    pseudobulk_qcs[pseudobulk_name].append(adata_xcs.obs)
    # aggregate across pseudobulks and save
    print("concat...")
    for p, x in pseudobulk_adatas.items():
        print(p)
        p_concat = ad.concat(x, axis=0)
        p_qc = pd.concat(pseudobulk_qcs[p], axis=0)
        # Save QC
        p_qc.to_csv(f"{data_dir}/rna_qc_reports/{p}-pseudobulked_cell_QC_metrics.csv", sep="\t", index=False)
        # Save h5ad
        p_concat.write(f"{data_dir}/pseudobulked_rna/{p}.h5ad")
        # make pseudobulk
        counts_df_p = p_concat.var.copy()
        counts_df_p["counts"] = p_concat.X.sum(axis=0).A1
        counts_df_p["CPM"] = (counts_df_p["counts"] / counts_df_p["counts"].sum()) * 1e6
        counts_df_p["log10CPM"] = np.log10(counts_df_p["CPM"] + 1)
        counts_df_p.to_csv(f"{data_dir}/pseudobulked_rna/{p}-pseudobulked_counts.tsv", sep="\t")


def main():
    parser = argparse.ArgumentParser(description="separate RNA h5ads file by pseudorep")

    # Define the expected flags
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')
    parser.add_argument('-g', '--geneinfo', type=str, required=True, help='Gene info file path')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    process_h5ad(args.datadir, args.metadata, args.geneinfo)
    

if __name__ == "__main__":
    main()

