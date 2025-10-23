import argparse
import json
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(description="Rename files")

    # Define the expected flags
    parser.add_argument('-o', '--outdir', type=str, required=True, help='Output file directory location')

    # Parse the arguments
    args = parser.parse_args()

    # Process file
    for x in os.listdir(f"{args.outdir}/jsons"):
        x_name = x[:-5]
        with open(f"{args.outdir}/jsons/{x}", 'r') as f:
            data = json.load(f)
        analysis_set = data["analysis_set"]
        counts_matrix = data["counts_matrix"]
        fragments = data["fragments"]
        assert(analysis_set == x_name)
        assert(os.path.exists(f"{args.outdir}/json_downloads/{counts_matrix}.h5ad"))
        assert(os.path.exists(f"{args.outdir}/json_downloads/{fragments}.bed.gz"))

        subprocess.run(["mv", f"{args.outdir}/json_downloads/{counts_matrix}.h5ad", f"{args.outdir}/raw_rna/{analysis_set}.h5ad"])
        subprocess.run(["mv", f"{args.outdir}/json_downloads/{fragments}.bed.gz", f"{args.outdir}/raw_fragments/{analysis_set}.bed.gz"])


if __name__ == "__main__":
    main()
