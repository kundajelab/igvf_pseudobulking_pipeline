import argparse
import os
import re

import pandas as pd


def remove_acgt_sequence(s):
    return re.sub(r'[ACGT]{16}', '', s)

def get_acgt_sequence(s):
    return re.search(r'[ACGT]{16}', s).group(0)

def decipher_lanes(annotations_file_loc, barcode_column, frags_dir, new_annotations_loc):
    # load fragment barcodes
    accession_barcodes = dict()
    for x in os.listdir(frags_dir):
        print(f"loading barcodes from {x.split('.')[0]}")
        df = pd.read_csv(f"{frags_dir}/{x}", sep="\t", compression="gzip", names=["chro", "start", "end", "barcode", "reads"], nrows=50000)
        df_barcodes = set([z.split("_")[0] for z in df["barcode"]])
        accession_barcodes[x.split(".")[0]] = df_barcodes
    # load annotations
    annotations_df = pd.read_csv(annotations_file_loc, sep="\t")
    annotations_df["barcode_just_dna"] = [get_acgt_sequence(x) for x in annotations_df[barcode_column]]
    annotations_df["barcode_without_dna"] = [remove_acgt_sequence(x) for x in annotations_df[barcode_column]]
    barcode_lanes = sorted(annotations_df["barcode_without_dna"].unique())
    # decipher
    assert len(barcode_lanes) == len(accession_barcodes), f"number of lane identifiers ({len(barcode_lanes)}) does not match number of accessions ({len(accession_barcodes)})"
    if len(barcode_lanes) == 1:
        print(f"trivial deciphering {barcode_lanes[0]} to {list(accession_barcodes.keys())[0]}")
        deciphering = {barcode_lanes[0]: list(accession_barcodes.keys())[0]}
    else:
        deciphering = dict()
        mapped_to = []
        for x in barcode_lanes:
            print(f"--- {x} ---")
            deciphering_x = []
            annotations_df_x = annotations_df[annotations_df["barcode_without_dna"] == x]
            barcodes_x = annotations_df_x["barcode_just_dna"]
            for a, b in accession_barcodes.items():
                b_match = sum(z in b for z in barcodes_x)
                print(f"* {a} {b_match}")
                deciphering_x.append((b_match, a))
            deciphering_x = sorted(deciphering_x, reverse=True)
            assert deciphering_x[0][0] > 2.5*deciphering_x[1][0], f"could not decipher lanetag {x}"
            assert deciphering_x[0][1] not in mapped_to, f"ERROR: {deciphering[0][1]} MAPPED TO MULTIPLE LANES"
            deciphering[x] = deciphering_x[0][1]
            print(f"{x} was mapped to {deciphering[x]}!")
    # create new annotations
    annotations_df["analysis_accession"] = [deciphering[x] for x in annotations_df["barcode_without_dna"]]
    annotations_df["barcode"] = annotations_df["barcode_just_dna"]
    annotations_df = annotations_df.drop(columns=["barcode_just_dna", "barcode_without_dna"])
    print(annotations_df)
    annotations_df.to_csv(new_annotations_loc, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="decipher which lane matches which lane code")

    # Define the expected flags
    parser.add_argument('-a', '--annotations', type=str, required=True, help='Annotations file')
    parser.add_argument('-bc', '--barcode_column', type=str, required=True, help='Barcode column')
    parser.add_argument('-fd', '--fragsdir', type=str, required=True, help='Fragments directory')
    parser.add_argument('-o', '--output', type=str, required=True, help='New annotations file save loc')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    decipher_lanes(args.annotations, args.barcode_column, args.fragsdir, args.output)
    

if __name__ == "__main__":
    main()

