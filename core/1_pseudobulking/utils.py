from numba import njit
import numpy as np
import pandas as pd


def load_metadata(metadata_loc):
    df = pd.read_csv(metadata_loc, sep="\t")
    # REQUIRED COLUMNS - barcode_sample, cell_name, cell_description, CL_id, CL_term_name, subsample, analysis_set_accession
    assert "barcode_sample" in df.columns, "metadata must contain 'barcode_sample' column"
    assert "cell_name" in df.columns, "metadata must contain 'cell_name' column"
    assert "cell_description" in df.columns, "metadata must contain 'cell_description' column"
    assert "CL_id" in df.columns, "metadata must contain 'CL_id' column"
    assert "CL_term_name" in df.columns, "metadata must contain 'CL_term_name' column"
    assert "subsample" in df.columns, "metadata must contain 'subsample' column"
    assert "analysis_set_accession" in df.columns, "metadata must contain 'analysis_set_accession' column"
    # COLUMN STRUCTURE REQUIREMENTS
    # TODO: SUBSAMPLE CANNOT CONTAIN HYPHENS OR CHARACTERS THAT WOULD BREAK FILE NAMING

    # RETURN
    return df


def map_cell_names_to_annotations(metadata_df, data_dir):
    # Cell name to annotation mapping
    cell_name_to_annotation_dict = {x: f"annotation_{i}" for i, x in enumerate(sorted(metadata_df["cell_name"].unique()))}
    # Map cell names to annotations in metadata_df
    metadata_df["annotation"] = metadata_df["cell_name"].map(cell_name_to_annotation_dict)
    # Save mapping to file
    cell_name_to_annotation_df = pd.DataFrame()
    cell_name_to_annotation_df["cell_name"] = list(cell_name_to_annotation_dict.keys())
    cell_name_to_annotation_df["annotation"] = list(cell_name_to_annotation_dict.values())
    cell_name_to_annotation_df.to_csv(f"{data_dir}/cell_name_to_annotation_mapping.tsv", sep="\t", index=False)
    # Return
    return metadata_df


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
