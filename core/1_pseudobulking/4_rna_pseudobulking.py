import argparse
from collections import defaultdict
import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

def load_metadata(metadata_loc):
    df = pd.read_csv(metadata_loc, sep="\t")
    assert "analysis_accession" in df.columns, "metadata must contain 'analysis_accession' column"
    assert "barcode" in df.columns, "metadata must contain 'barcode' column"
    assert "annotation" in df.columns, "metadata must contain 'annotation' column"
    assert "subsample" in df.columns, "metadata must contain 'subsample' column"
    assert all(["-" not in x for x in df["subsample"].unique().tolist()]), "'subsample' column cannot contain '-' chracter"
    for c in df.columns:
        if c.startswith("annotation"):
            annotations_c = df[c].unique().tolist()
            assert all(["-" not in x for x in annotations_c]), f"'{c}' column cannot contain '-' character"
    return df

def process_h5ad(data_dir, metadata_loc):
    # Load metadata
    metadata_df = load_metadata(metadata_loc)
    print(metadata_df)
    # Iterate through h5ads
    pseudobulk_adatas = defaultdict(list)
    qc_df_list = []
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    gene_ref = pd.read_csv(f"{repo_root}/chr_info_data/GRCh38_GRCm39_gene_id_to_name.csv", index_col=0) ## convert the path
    for x in os.listdir(f"{data_dir}/raw_rna"):
        # Load h5ads
        x_name = x.split(".")[0]
        metadata_df_x = metadata_df[metadata_df["analysis_accession"] == x_name]
        adata = ad.read_h5ad(f"{data_dir}/raw_rna/{x}")
        print(x_name)
        print(adata)
        # Compute QC for each file (analysis_accession)
        adata.var['gene_symbol'] = adata.var.index.map(gene_ref['gene_name'])
        adata.var['mt']          = adata.var.index.map(gene_ref['mt']).fillna(False).astype(bool)
        adata.var['ribo']        = adata.var.index.map(gene_ref['ribo']).fillna(False).astype(bool)
        sc.pp.calculate_qc_metrics(
            adata, 
            qc_vars=['mt', 'ribo'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        adata.obs['barcode'] = adata.obs.index
        qc_df_sub = adata.obs[['barcode',"total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]].copy()
        qc_df_sub.columns = ['barcode', 'read_count', 'gene_count', 'pct_mito', 'pct_ribo']
        qc_df_sub['read_count'] = qc_df_sub['read_count'].astype(int)
        qc_df_sub.to_csv(f"{data_dir}/rna_qc_reports/{x_name}_scRNA_total_QC_metrics.csv", index=False) # save all cells QC results for each analysis accession
        qc_df_sub_filtered = qc_df_sub[qc_df_sub['barcode'].isin(metadata_df_x['barcode'])]
        qc_df_list.append(qc_df_sub_filtered) # only keep barcodes in metadata (annotation file) to generate concat QC files across all analysis accessions
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
                    adata_xcs.obs["analysis_accession"] = x_name
                    pseudobulk_adatas[pseudobulk_name].append(adata_xcs)
    qc_df = pd.concat(qc_df_list, axis=0)
    # aggregate across pseudobulks and save
    print("concat...")
    for p, x in pseudobulk_adatas.items():
        print(p)
        p_concat = ad.concat(x, axis=0)
        # save QC file matching pseudobulked barcodes
        p_qc_df = qc_df[qc_df['barcode'].isin(p_concat.obs_names)]
        p_qc_df.to_csv(f"{data_dir}/rna_qc_reports/{p}_pseudobulked_QC_metrics.csv", index=False)
        # save h5ad
        p_concat.write(f"{data_dir}/pseudobulked_rna/{p}.h5ad")
        # make pseudobulk
        counts_df_p = p_concat.var.copy()
        counts_df_p["counts"] = p_concat.X.sum(axis=0).A1
        counts_df_p.to_csv(f"{data_dir}/pseudobulked_rna/{p}-pseudobulked_counts.tsv", sep="\t")



def main():
    parser = argparse.ArgumentParser(description="separate fragments file by pseudorep")

    # Define the expected flags
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    process_h5ad(args.datadir, args.metadata)
    

if __name__ == "__main__":
    main()

