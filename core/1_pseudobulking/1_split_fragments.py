import argparse
import pathlib
import random

import gzip
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


def process_fragments_file(fragments_file_name, data_dir, metadata_loc, chr_order_loc):
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
    chr_order_df = pd.read_csv(chr_order_loc, sep="\t", names=["chr", "size"])
    allowed_chrs = set(chr_order_df["chr"].unique().tolist())
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
    num_lines_notfound = 0
    fragments_file_loc = f"{data_dir}/raw_fragments/{fragments_file_name}.bed.gz"
    with gzip.open(fragments_file_loc, 'rt') as f:
        for line in f:
            # Parse line
            num_lines += 1
            chro, start, end, barcode, reads = tuple(line.strip().split("\t"))
            start, end = int(start), int(end)
            # TODO: SHIFTING
            start += 0
            end += 0
            # Remove nonstandard chromosomes
            if chro not in allowed_chrs:
                num_lines_notfound += 1
                continue
            # Put line into new pseudobulks
            line_pseudobulks = barcodes_to_pseudobulks.get(barcode, None)
            if line_pseudobulks is None:
                num_lines_notfound += 1
                continue
            else:
                num_lines_found += 1
                for pseudobulk in line_pseudobulks:
                    # Write fragment to fragments
                    output_file_handles["fragments"][pseudobulk].write(line)
                    # Write insertions to pseudorepT
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode}\t{reads}\n")
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode}\t{reads}\n")
                    # Write start insertion to pseudorep1/2
                    if random.randint(0, 1) == 0:
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode}\t{reads}\n")
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode}\t{reads}\n")
                    # Write end insertion to pseudorep1/2
                    if random.randint(0, 1) == 0:
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode}\t{reads}\n")
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode}\t{reads}\n")
        # close outfiles
        for f in output_file_handles["pseudorep1"].values():
            f.close()
        for f in output_file_handles["pseudorep2"].values():
            f.close()
        for f in output_file_handles["pseudorepT"].values():
            f.close()
        for f in output_file_handles["fragments"].values():
            f.close()
        # results
        print(f"{fragments_file_name}: {num_lines} fragments --> {num_lines_found} kept/{num_lines_notfound} filtered")


def main():
    parser = argparse.ArgumentParser(description="separate fragments file by pseudorep")

    # Define the expected flags
    parser.add_argument('-f', '--fragments', type=str, required=True, help='Input fragments file path')
    parser.add_argument('-d', '--datadir', type=str, required=True, help='Data directory')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Input annotations metadata file path')
    parser.add_argument('-c', '--chr_order', type=str, required=True, help='Chromosomes order/size file path')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    process_fragments_file(args.fragments, args.datadir, args.metadata, args.chr_order)
    

if __name__ == "__main__":
    main()

