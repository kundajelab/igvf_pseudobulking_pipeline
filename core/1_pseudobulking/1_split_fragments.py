import argparse
from collections import defaultdict
import random
import time

import gzip
import numpy as np
import pandas as pd

from utils import load_metadata, get_pseudobulk_name, load_tss_locs, check_tss_overlap


##################
# MAIN FUNCTIONS #
##################
def process_fragments_file(fragments_file_name, data_dir, metadata_loc, at_annotation_level, chr_sizes_loc, tss_locs_loc):
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
        row_pseudobulks = [get_pseudobulk_name(x, row[x], row_sample, at_annotation_level) for x in annotation_columns]
        barcodes_to_pseudobulks[row_barcode] = row_pseudobulks
        pseudobulks.update(row_pseudobulks)
    # Get allowed chromosomes
    chr_sizes_df = pd.read_csv(chr_sizes_loc, sep="\t", names=["chr", "size"])
    allowed_chrs = set(chr_sizes_df["chr"].unique().tolist())
    # Get TSS locations/data
    tss_by_chr_np = load_tss_locs(tss_locs_loc)
    TSS_half_window = 2000
    TSS_half_smooth_window = 5
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
    fragments_file_loc = f"{data_dir}/raw_fragments/{fragments_file_name}.bed.gz"
    cell_qc = defaultdict(lambda: {"annotated": False, "num_unique_frags": 0, "num_reads": 0, "num_dup_reads": 0, "mono_nucleosomal_frags": 0, "nucleosome_free_frags": 0, "tss_insertions": np.zeros((2*TSS_half_window+1, ))})
    num_lines = 0
    start_time = time.time()
    with gzip.open(fragments_file_loc, 'rt') as f:
        for line in f:
            # if num_lines >= 5e6:
            #     break
            # Parse line
            num_lines += 1
            chro, start, end, barcode, reads = tuple(line.strip().split("\t"))
            start, end, reads = int(start), int(end), int(reads)
            # Shifting (+4/-4)
            start_shifted = start + 4
            end_shifted = end - 4
            # Cell QC
            barcode_qc = cell_qc[barcode]
            barcode_qc["num_unique_frags"] += 1
            barcode_qc["num_reads"] += reads
            barcode_qc["num_dup_reads"] += reads-1
            # nucleosomal enrichment
            barcode_qc["nucleosome_free_frags"] += int((end_shifted - start_shifted) < 148)
            barcode_qc["mono_nucleosomal_frags"] += int((148 <= (end_shifted - start_shifted) < 295))
            # tss enrichment
            if (chro in tss_by_chr_np) and (chro not in ["chrM"]):
               start_tss_pos = check_tss_overlap(start_shifted, tss_by_chr_np[chro][0], tss_by_chr_np[chro][1])
               if start_tss_pos is not None:
                   for p in start_tss_pos:
                       # p is between -L and L
                       barcode_qc["tss_insertions"][p+TSS_half_window] += 1
               end_tss_pos = check_tss_overlap(end_shifted-1, tss_by_chr_np[chro][0], tss_by_chr_np[chro][1])
               if end_tss_pos is not None:
                   for p in end_tss_pos:
                       # p is between -L and L
                       barcode_qc["tss_insertions"][p+TSS_half_window] += 1
            # Remove nonstandard chromosomes
            if chro not in allowed_chrs:
                continue
            # Put line into new pseudobulks
            line_pseudobulks = barcodes_to_pseudobulks.get(barcode, None)
            if line_pseudobulks is None:
                continue
            else:
                barcode_qc["annotated"] = True
                for pseudobulk in line_pseudobulks:
                    # Write fragment to fragments
                    output_file_handles["fragments"][pseudobulk].write(line)
                    # Get start/end insertions
                    start_insertion_line = f"{chro}\t{start_shifted}\t{start_shifted+1}\t{barcode}\t{1}\n"
                    end_insertion_line = f"{chro}\t{end_shifted-1}\t{end_shifted}\t{barcode}\t{1}\n"
                    # Write insertions to pseudorepT
                    output_file_handles["pseudorepT"][pseudobulk].write(start_insertion_line)
                    output_file_handles["pseudorepT"][pseudobulk].write(end_insertion_line)
                    # Write insertions to pseudorep1/2
                    if random.random() < 0.5:
                        output_file_handles["pseudorep1"][pseudobulk].write(start_insertion_line)
                        output_file_handles["pseudorep1"][pseudobulk].write(end_insertion_line)
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(start_insertion_line)
                        output_file_handles["pseudorep2"][pseudobulk].write(end_insertion_line)
    end_time = time.time()
    print(f"({fragments_file_name}) Time taken: {end_time-start_time:.3f} seconds for {num_lines} lines ({num_lines/(end_time-start_time):.3f} lines/second)")
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
    tss_rows = []
    for barcode, barcode_qc in cell_qc.items():
        qc_row = dict()
        # Summary values
        qc_row["analysis_accession"] = fragments_file_name
        qc_row["barcode"] = barcode
        qc_row["annotated"] = barcode_qc["annotated"]
        qc_row["num_frags"] = barcode_qc["num_unique_frags"]
        qc_row["pct_duplicated_reads"] = (barcode_qc["num_dup_reads"]/barcode_qc["num_reads"]) * 100
        qc_row["nucleosomal_signal"] = (1+barcode_qc["mono_nucleosomal_frags"])/(1+barcode_qc["nucleosome_free_frags"])
        tss_insertions = barcode_qc["tss_insertions"]
        tss_insertions_flank_mean = (np.sum(tss_insertions[:100]) + np.sum(tss_insertions[-100:]))/200
        tss_insertions_center = np.mean(tss_insertions[TSS_half_window-TSS_half_smooth_window:TSS_half_window+TSS_half_smooth_window+1])
        qc_row["tss_enrichment"] = tss_insertions_center/(tss_insertions_flank_mean + 0.1) # add 0.1 like snapatac2 to avoid division by zero
        # Raw values
        qc_row["raw-num_reads"] = barcode_qc["num_reads"]
        qc_row["raw-num_dup_reads"] = barcode_qc["num_dup_reads"]
        qc_row["raw-mono_nucleosomal_frags"] = barcode_qc["mono_nucleosomal_frags"]
        qc_row["raw-nucleosome_free_frags"] = barcode_qc["nucleosome_free_frags"]
        # Append row
        qc_rows.append(qc_row)
        tss_rows.append(barcode_qc["tss_insertions"])
    df = pd.DataFrame(qc_rows)
    df.to_csv(f"{data_dir}/atac_qc_reports/{fragments_file_name}.tsv", sep="\t", index=False)
    tss_matrix = np.vstack(tss_rows)
    np.save(f"{data_dir}/atac_qc_reports/{fragments_file_name}_tss_matrix.npy", tss_matrix)


def main():
    parser = argparse.ArgumentParser(description="separate fragments file by pseudorep")

    # Define the expected flags
    parser.add_argument('-f', '--fragments', type=str, required=True, help='Input fragments file path')
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')
    parser.add_argument('-a', '--at_annotation_level', type=str, required=True, help='At annotation level flag')
    parser.add_argument('-c', '--chr_sizes', type=str, required=True, help='Chromosomes size/order file path')
    parser.add_argument('-t', '--tss_locs', type=str, required=True, help="TSS locations file")

    # Parse the arguments
    args = parser.parse_args()

    if args.at_annotation_level == "True":
        args.at_annotation_level = True
    elif args.at_annotation_level == "False":
        args.at_annotation_level = False
    else:
        raise ValueError("at_annotation_level must be 'True' or 'False'")

    # Process file
    process_fragments_file(args.fragments, args.datadir, args.metadata, args.at_annotation_level, args.chr_sizes, args.tss_locs)
    

if __name__ == "__main__":
    main()

