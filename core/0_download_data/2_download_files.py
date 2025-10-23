import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description="Download IGVF files")

    # Define the expected flags
    parser.add_argument('-k', '--access_key', type=str, required=True, help='IGVF download access key')
    parser.add_argument('-s', '--secret_key', type=str, required=True, help='IGVF download secret key')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file directory location')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    for x in os.listdir(f"{args.outdir}/jsons"):
        subprocess.run(["bash", f"{current_dir}/download_files.sh", "--access_key", args.access_key, "--secret_key", args.secret_key, "--igvffs_in", f"{args.outdir}/jsons/{x}", "--outdir", f"{args.outdir}/json_downloads"])


if __name__ == "__main__":
    main()
