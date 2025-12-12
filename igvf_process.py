import argparse
import os
import subprocess

import pandas as pd


def get_species_files(species):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    species_dir = os.path.join(scriptdir, "genome_data", species)
    if not os.path.exists(species_dir):
        raise ValueError(f"Invalid species {species}. Must be either human or mouse.")
    chr_sizes = f"{species_dir}/chr_sizes.tsv"
    blacklist_file = f"{species_dir}/blacklist.bed"
    tss_file = f"{species_dir}/tss.tsv"
    gene_info_file = f"{species_dir}/gene_info.csv"
    return chr_sizes, blacklist_file, tss_file, gene_info_file


def main():
    parser = argparse.ArgumentParser(description="IGVF single cell pseudobulking pipeline")

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # download
    download_parser = subparsers.add_parser("download", help="Download uniform pipeline outputs by uniform analysis accession")
    download_parser.add_argument("-al", required=True, help="Accession list (comma separated)")
    download_parser.add_argument("-w", required=True, help="Workspace")
    download_parser.add_argument("-ak", required=True, help="Access key")
    download_parser.add_argument("-sk", required=True, help="Secret key")

    # pseudobulk
    pseudobulk_parser = subparsers.add_parser("pseudobulk", help="Pseudobulk data given uniform annotation file")
    pseudobulk_parser.add_argument("-w", required=True, help="Workspace")
    pseudobulk_parser.add_argument("-a", required=True, help="Annotation file")
    pseudobulk_parser.add_argument("-s", required=True, help="Species")
    pseudobulk_parser.add_argument("-c", required=True, help="Number of CPUs")
    pseudobulk_parser.add_argument("--at_annotation_level", action="store_true", help="If set, pseudobulk at annotation level instead of annotation x subsample level")

    # full
    full_parser = subparsers.add_parser("full", help="Download + pseudobulk data given uniform annotation file")
    full_parser.add_argument("-a", required=True, help="Annotation file")
    full_parser.add_argument("-w", required=True, help="Workspace")
    full_parser.add_argument("-s", required=True, help="Species")
    full_parser.add_argument("-ak", required=True, help="Access key")
    full_parser.add_argument("-sk", required=True, help="Secret key")
    full_parser.add_argument("-c", required=True, help="Number of CPUs")
    full_parser.add_argument("--at_annotation_level", action="store_true", help="If set, pseudobulk at annotation level instead of annotation x subsample level")

    # decipher
    decipher_parser = subparsers.add_parser("decipher", help="Decipher lane to accession map given non-uniform annotation file")
    decipher_parser.add_argument("-a", required=True, help="Annotation file")
    decipher_parser.add_argument("-blc", required=True, help="barcode_lane column (Column name in annotation file containing barcode_lane values)")
    decipher_parser.add_argument("-w", required=True, help="Workspace")

    # Parse the arguments
    args = parser.parse_args()

    # Check if a legal command was given
    if not args.command:
        parser.print_help()
        exit(1)

    # Run pipeline based on command
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    if args.command == "download":
        download_path = os.path.join(scriptdir, "core", "0_download_data", "PIPELINE.sh")
        download_command = ["bash", download_path, args.w, args.al, args.ak, args.sk]
        subprocess.run(download_command)
    
    elif args.command == "pseudobulk":
        chr_sizes, blacklist_file, tss_file, gene_info_file = get_species_files(args.s)

        pseudobulk_path = os.path.join(scriptdir, "core", "1_pseudobulking", "PIPELINE.sh")
        pseudobulk_command = ["bash", pseudobulk_path, args.w, args.a, str(args.at_annotation_level), chr_sizes, blacklist_file, tss_file, gene_info_file, args.c]
        subprocess.run(pseudobulk_command)
    
    elif args.command == "full":
        annotation_df = pd.read_csv(args.a, sep="\t")
        annotation_list = ",".join(annotation_df["analysis_accession"].unique().tolist())

        download_path = os.path.join(scriptdir, "core", "0_download_data", "PIPELINE.sh")
        download_command = ["bash", download_path, args.w, annotation_list, args.ak, args.sk]
        subprocess.run(download_command)

        chr_sizes, blacklist_file, tss_file, gene_info_file = get_species_files(args.s)

        pseudobulk_path = os.path.join(scriptdir, "core", "1_pseudobulking", "PIPELINE.sh")
        pseudobulk_command = ["bash", pseudobulk_path, args.w, args.a, str(args.at_annotation_level), chr_sizes, blacklist_file, tss_file, gene_info_file, args.c]
        subprocess.run(pseudobulk_command)
    
    elif args.command == "decipher":
        decipher_path = os.path.join(scriptdir, "core", "decipher", "PIPELINE.sh")
        decipher_command = ["bash", decipher_path, args.a, args.blc, args.w]
        subprocess.run(decipher_command)


if __name__ == "__main__":
    main()
