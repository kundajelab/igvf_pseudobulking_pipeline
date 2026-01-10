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


def get_pseudobulk_name(annotation_type, annotation, subsample, at_annotation_level):
    if at_annotation_level:
        pseudobulk_name = f"{annotation_type}-{annotation}"
    else:
        pseudobulk_name = f"{annotation_type}-{annotation}-{subsample}"
    return pseudobulk_name


def load_tss_locs(tss_locs_loc):
    tss_locs_df = pd.read_csv(tss_locs_loc, sep="\t")
    tss_by_chr_np = dict()
    for chro in tss_locs_df["chro"].unique():
        chro_tss_df = tss_locs_df[tss_locs_df["chro"] == chro]
        tss_positions = chro_tss_df["TSS"].values - 1 # GTF is 1-based but fragments are 0-based
        strand_signs = np.array([1 if s == "+" else -1 for s in chro_tss_df["strand"].values], dtype=int)
        df = pd.DataFrame()
        df["TSS"] = tss_positions
        df["strand"] = strand_signs
        df.sort_values(by=["TSS", "strand"], inplace=True)
        df = df.drop_duplicates()
        tss_by_chr_np[chro] = (df["TSS"].to_numpy(), df["strand"].to_numpy())
    return tss_by_chr_np


@njit
def check_tss_overlap(position, tss_vec, strand_vec):
    L = 2000 # half window size
    tss_left_idx = np.searchsorted(tss_vec, position-L, side='left')
    tss_right_idx = np.searchsorted(tss_vec, position+L, side='right')-1
    if tss_left_idx > tss_right_idx:
        return None
    return [(position - tss_vec[i]) * strand_vec[i] for i in range(tss_left_idx, tss_right_idx+1)]
