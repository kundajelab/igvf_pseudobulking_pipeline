import argparse
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(description="Download IGVF files")

    # Define the expected flags
    parser.add_argument('-a', '--accessions', type=str, required=True, help='IGVF accessions')
    parser.add_argument('-k', '--access_key', type=str, required=True, help='IGVF download access key')
    parser.add_argument('-s', '--secret_key', type=str, required=True, help='IGVF download secret key')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file directory location')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    accessions = args.accessions.split(",")
    for a in accessions:
        subprocess.run(["python", f"{current_dir}/get_igvf_download_accessions.py", "-a", a, "--access_key", args.access_key, "--secret_key", args.secret_key, "-o", f"{args.outdir}/jsons/{a}.json"])


if __name__ == "__main__":
    main()
