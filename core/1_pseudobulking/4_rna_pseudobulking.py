import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd


def load_metadata(metadata_loc):
    df = pd.read_csv(metadata_loc, sep="\t")
    assert "barcode" in df.columns, "metadata must contain 'barcode' column"
    assert "analysis_accession" in df.columns, "metadata must contain 'analysis_accession' column"
    return df


def process_h5ad(data_dir, metadata_loc):
    # load metadata
    df = load_metadata(metadata_loc)
    # get all pseudobulks
    pseudobulks = []
    for x in df.columns:
        if x in ["barcode", "analysis_accession"]:
            continue
        pseudobulks_x = df[x].unique().tolist()
        assert all([z not in pseudobulks for z in pseudobulks_x]), "all pseudobulks must have unique names"
        pseudobulks += pseudobulks_x
    assert all(["-" not in z for z in pseudobulks]), "no '-' in pseudobulk names"
    # iterate through h5ads
    pseudobulk_adatas = {x: [] for x in pseudobulks}
    for x in os.listdir(f"{data_dir}/raw_rna"):
        x_name = x.split(".")[0]
        df_x = df[df["analysis_accession"] == x_name]
        adata = ad.read_h5ad(f"{data_dir}/raw_rna/{x}")
        for c in df_x.columns:
            if c in ["barcode", "analysis_accession"]:
                continue
            df_x_c_revdict = {row["barcode"]: row[c] for _, row in df_x.iterrows()}
            adata.obs[c] = [df_x_c_revdict.get(b.split("_")[0], "NOTINDICT") for b in adata.obs_names]
        for c in df_x.columns:
            if c in ["barcode", "analysis_accession"]:
                continue
            for p in df_x[c].unique():
                adata_p = adata[adata.obs[c] == p, :].copy()
                adata_p.obs["uniform_analysis_accession"] = x_name
                adata_p.obs = adata_p.obs[["uniform_analysis_accession"]]
                pseudobulk_adatas[p].append(adata_p)
    # aggregate across pseudobulks and save
    for p, x in pseudobulk_adatas.items():
        p_concat = ad.concat(x, axis=0)
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

