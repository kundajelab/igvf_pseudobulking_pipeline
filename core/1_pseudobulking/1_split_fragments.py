import argparse
import pathlib
import random

import gzip
import pandas as pd


def load_metadata(metadata_loc):
    df = pd.read_csv(metadata_loc, sep="\t")
    assert "barcode" in df.columns, "metadata must contain 'barcode' column"
    assert "analysis_accession" in df.columns, "metadata must contain 'analysis_accession' column"
    return df


def process_fragments_file(fragments_file_name, data_dir, metadata_loc, chr_order_loc):
    assert "-" not in fragments_file_name, "no '-' in fragments file name"
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
    # set up barcodes to pseudobulks
    barcodes_to_pseudobulks = dict()
    analysis_accessions = df["analysis_accession"]
    for x in analysis_accessions:
        barcodes_to_pseudobulks[x] = dict()
    for _, row in df.iterrows():
        row_pseudobulks = [y for x, y in row.items() if x not in ["barcode", "analysis_accession"]]
        barcodes_to_pseudobulks[row['analysis_accession']][row['barcode']] = row_pseudobulks
    barcodes = {x: set(barcodes_to_pseudobulks[x].keys()) for x in barcodes_to_pseudobulks}
    # get allowed chromosomes
    chr_order_df = pd.read_csv(chr_order_loc, sep="\t", names=["chr", "size"])
    allowed_chrs = chr_order_df["chr"].unique().tolist()
    # set up outfiles
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
    # iterate through fragments file
    num_lines = 0
    num_lines_found = 0
    num_lines_notfound = 0
    fragments_file_loc = f"{data_dir}/raw_fragments/{fragments_file_name}.bed.gz"
    with gzip.open(fragments_file_loc, 'rt') as f:
        for line in f:
            # Parse line
            num_lines += 1
            chro, start, end, barcode_accession, reads = tuple(line.strip().split("\t"))
            barcode, accession = tuple(barcode_accession.split("_"))
            start, end = int(start), int(end)
            if barcode in barcodes[fragments_file_name] and chro in allowed_chrs:
                num_lines_found += 1
                # Look through all pseudobulks
                for pseudobulk in barcodes_to_pseudobulks[fragments_file_name][barcode]:
                    # Write fragment to fragments
                    output_file_handles["fragments"][pseudobulk].write(line)
                    # Write insertions to pseudorepT
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode_accession}\t{reads}\n")
                    output_file_handles["pseudorepT"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode_accession}\t{reads}\n")
                for pseudobulk in barcodes_to_pseudobulks[fragments_file_name][barcode]:
                    # Write start insertion
                    if random.randint(0, 1) == 0:
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode_accession}\t{reads}\n")
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{start}\t{start+1}\t{barcode_accession}\t{reads}\n")
                    # Write end insertion
                    if random.randint(0, 1) == 0:
                        output_file_handles["pseudorep1"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode_accession}\t{reads}\n")
                    else:
                        output_file_handles["pseudorep2"][pseudobulk].write(f"{chro}\t{end-1}\t{end}\t{barcode_accession}\t{reads}\n")
            else:
                # Did not find line
                num_lines_notfound += 1
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

