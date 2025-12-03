import argparse
import random
import time

import gzip
from numba import njit
import numpy as np
import pandas as pd


####################
# HELPER FUNCTIONS #
####################
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
        tss_by_chr[row["chro"]][row["transcript"]] = row["TSS"]
    tss_by_chr_np = {x: np.array(list(tss_by_chr[x].values()), dtype=int) for x in tss_by_chr}
    return tss_by_chr_np


@njit
def check_tss_overlap(position, tss_vec):
    for i in range(len(tss_vec)):
        distance = tss_vec[i] - position
        if ((-1000 <= distance) and (distance <= 1000)):
            return distance+1000 # return first TSS distance found
    return None


##################
# MAIN FUNCTIONS #
##################
def process_fragments_file(fragments_file_name, data_dir, metadata_loc, chr_sizes_loc, tss_locs_loc):
    # Load metadata
    metadata_df = load_metadata(metadata_loc)
    # Subset metadata to current analysis accession (fragment file name)
    metadata_df = metadata_df[metadata_df["analysis_accession"] == fragments_file_name]
    # Get annotation columns
    annotation_columns = [x for x in metadata_df.columns if x.startswith("annotation")]
    # Compute barcodes --> annotation mapping
    pseudobulks = set() # {annotation_type}-{annotation}-{subsample}
    barcodes_to_pseudobulks = dict()
    for _, row in metadata_df.iterrows():
        row_barcode = row["barcode"]
        row_sample = row["subsample"]
        row_pseudobulks = [f"{x}-{row[x]}-{row_sample}" for x in annotation_columns]
        barcodes_to_pseudobulks[row_barcode] = row_pseudobulks
        pseudobulks.update(row_pseudobulks)
    # Get allowed chromosomes
    chr_sizes_df = pd.read_csv(chr_sizes_loc, sep="\t", names=["chr", "size"])
    allowed_chrs = set(chr_sizes_df["chr"].unique().tolist())
    # Get TSS locations
    tss_by_chr_np = load_tss_locs(tss_locs_loc)
    # Set up outfiles
    output_file_handles = dict()
    output_file_handles["pseudorep1"] = dict()
    output_file_handles["pseudorep2"] = dict()
    output_file_handles["pseudorepT"] = dict()
    output_file_handles["fragments"] = dict()
    for p in pseudobulks:
        output_file_handles["pseudorep1"][p] = open(f"{data_dir}/separated_pseudorep1/{p}-{fragments_file_name}.tsv", 'w')
        output_file_handles["pseudorep2"][p] = open(f"{data_dir}/separated_pseudorep2/{p}-{fragments_file_name}.tsv", 'w')
        output_file_handles["pseudorepT"][p] = open(f"{data_dir}/separated_pseudorepT/{p}-{fragments_file_name}.tsv", 'w')
        output_file_handles["fragments"][p] = open(f"{data_dir}/separated_fragments/{p}-{fragments_file_name}.tsv", 'w')
    # Iterate through fragments file
    num_lines = 0
    num_lines_found = 0
    fragments_file_loc = f"{data_dir}/raw_fragments/{fragments_file_name}.bed.gz"
    cell_qc = dict()
    start_time = time.time()
    with gzip.open(fragments_file_loc, 'rt') as f:
        for line in f:
            if (num_lines > 5000000):
                break
            # Parse line
            num_lines += 1
            chro, start, end, barcode, reads = tuple(line.strip().split("\t"))
            start, end, reads = int(start), int(end), int(reads)
            
            # Shifting (+4/-4)
            start_shifted = start + 4
            end_shifted = end - 4
            # Count fragment statistics
            if barcode not in cell_qc:
                cell_qc[barcode] = {"annotated": False, "num_frags": 0, "num_dup_reads": 0, "mono_nucleosomal_frags": 0, "nucleosome_free_frags": 0, "nonstandard_chr_frags": 0, "tss_insertions": np.zeros((2001, ))}
            # basic fragment counts
            cell_qc[barcode]["num_frags"] += reads
            cell_qc[barcode]["num_dup_reads"] += reads-1
            # nucleosomal enrichment
            frag_len = end-start
            if frag_len < 148:
                cell_qc[barcode]["nucleosome_free_frags"] += reads
            elif frag_len < 295:
                cell_qc[barcode]["mono_nucleosomal_frags"] += reads
            # tss enrichment
            if chro != "chrM" and chro in tss_by_chr_np:
                start_tss_pos = check_tss_overlap(start_shifted, tss_by_chr_np[chro])
                if start_tss_pos is not None:
                    cell_qc[barcode]["tss_insertions"][start_tss_pos] += 1
                end_tss_pos = check_tss_overlap(end_shifted, tss_by_chr_np[chro])
                if end_tss_pos is not None:
                    cell_qc[barcode]["tss_insertions"][end_tss_pos] += 1
            # Remove nonstandard chromosomes
            if chro not in allowed_chrs:
                cell_qc[barcode]["nonstandard_chr_frags"] += reads
                continue
            # Put line into new pseudobulks
            line_pseudobulks = barcodes_to_pseudobulks.get(barcode, None)
            if line_pseudobulks is None:
                continue
            else:
                cell_qc[barcode]["annotated"] = True
                num_lines_found += 1
                for pseudobulk in line_pseudobulks:
                    # Write fragment to fragments
                    output_file_handles["fragments"][pseudobulk].write(line)
                    # Write insertions to pseudorepT
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{start_shifted}\t{start_shifted+1}\t{barcode}\t{reads}\n")
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{end_shifted-1}\t{end_shifted}\t{barcode}\t{reads}\n")
                    # Write insertions to pseudorep1/2
                    if random.randint(0, 1) == 0:
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{start_shifted}\t{start_shifted+1}\t{barcode}\t{reads}\n")
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{end_shifted-1}\t{end_shifted}\t{barcode}\t{reads}\n")
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{start_shifted}\t{start_shifted+1}\t{barcode}\t{reads}\n")
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{end_shifted-1}\t{end_shifted}\t{barcode}\t{reads}\n")
            
    end_time = time.time()
    print(f"Time taken: {end_time-start_time:.3f} seconds")
    # Close outfiles
    for f in output_file_handles["pseudorep1"].values():
        f.close()
    for f in output_file_handles["pseudorep2"].values():
        f.close()
    for f in output_file_handles["pseudorepT"].values():
        f.close()
    for f in output_file_handles["fragments"].values():
        f.close()
    # QC results
    qc_rows = []
    for barcode, barcode_qc in cell_qc.items():
        qc_row = dict()
        qc_row["barcode"] = barcode
        qc_row["annotated"] = barcode_qc["annotated"]
        qc_row["num_frags"] = barcode_qc["num_frags"]
        qc_row["percent_duplicated_reads"] = barcode_qc["num_dup_reads"]/barcode_qc["num_frags"]
        qc_row["nucleosomal_signal"] = (1+barcode_qc["mono_nucleosomal_frags"])/(1+barcode_qc["nucleosome_free_frags"])
        tss_insertions = barcode_qc["tss_insertions"]
        tss_insertions_flank_mean = (np.sum(tss_insertions[:100]) + np.sum(tss_insertions[-100:]))/200
        tss_insertions *= (1/tss_insertions_flank_mean if tss_insertions_flank_mean > 0 else 0)
        qc_row["tss_enrichment_max"] = np.max(tss_insertions)
        qc_row["tss_enrichment_center"] = tss_insertions[1000]
        qc_row["tss_enrichment_mean"] = np.mean(tss_insertions[500:1501])
        qc_row["percent_nonstandard_chr_frags"] = barcode_qc["nonstandard_chr_frags"]/barcode_qc["num_frags"]
        qc_rows.append(qc_row)
    df = pd.DataFrame(qc_rows)
    df.to_csv(f"{data_dir}/atac_qc_reports/{fragments_file_name}.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description="separate fragments file by pseudorep")

    # Define the expected flags
    parser.add_argument('-f', '--fragments', type=str, required=True, help='Input fragments file path')
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')
    parser.add_argument('-c', '--chr_sizes', type=str, required=True, help='Chromosomes size/order file path')
    parser.add_argument('-t', '--tss_locs', type=str, required=True, help="TSS locations file")

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    process_fragments_file(args.fragments, args.datadir, args.metadata, args.chr_sizes, args.tss_locs)
    

if __name__ == "__main__":
    main()

