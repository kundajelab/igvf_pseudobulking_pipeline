from numba import njit
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


def load_tss_locs(tss_locs_loc):
    tss_locs_df = pd.read_csv(tss_locs_loc, sep="\t")
    tss_by_chr = {x: dict() for x in tss_locs_df["chro"].unique()}
    for _, row in tss_locs_df.iterrows():
        strand_sign = 1 if row["strand"] == "+" else -1
        tss_by_chr[row["chro"]][row["transcript"]] = (row["TSS"], strand_sign)
    tss_by_chr_np = {x: (np.array([v[0] for v in tss_by_chr[x].values()], dtype=int), np.array([v[1] for v in tss_by_chr[x].values()], dtype=int)) for x in tss_by_chr}
    return tss_by_chr_np


@njit
def check_tss_overlap(position, tss_vec, strand_vec):
    for i in range(len(tss_vec)):
        distance = strand_vec[i]*(tss_vec[i] - position)
        if ((-1000 <= distance) and (distance <= 1000)):
            return distance+1000 # return first TSS distance found
    return None