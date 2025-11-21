import argparse
from collections import defaultdict
import os

import anndata as ad
import numpy as np
import pandas as pd


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
    for x in os.listdir(f"{data_dir}/raw_rna"):
        # Load h5ads
        x_name = x.split(".")[0]
        metadata_df_x = metadata_df[metadata_df["analysis_accession"] == x_name]
        adata = ad.read_h5ad(f"{data_dir}/raw_rna/{x}")
        print(x_name)
        print(adata)
        # TODO: Compute QC
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
    # aggregate across pseudobulks and save
    print("concat...")
    for p, x in pseudobulk_adatas.items():
        print(p)
        p_concat = ad.concat(x, axis=0)
        # TODO: QC
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

