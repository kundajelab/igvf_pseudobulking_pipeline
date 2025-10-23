###
# WRITTEN BY VIVIAN HECHT
###

import argparse
import json
from pathlib import Path
import requests
from typing import Dict


def get_igvf_ids_for_download(
        analysis_accession: str,
        access_key: str,
        secret_key: str,
        base_url: str,
        outfile: Path) -> None:
    """Given an analysis IGVF ID, get the ids of the ATAC fragment file and
    the RNA-seq h5ad file and save to a json

    Args:
    analysis_accession: ID from IGVF uniform analysis pipeline
    access_key: IGVF access key
    secret_key: IGVF secret key
    base_url: URL of IGVF portal
    outfile: where to save output IDs
    """
    headers = {'accept': 'application/json'}
    igvf_auth = (access_key, secret_key)
    json_url = f"{base_url}/{analysis_accession}/?format=json"
    json_results = requests.get(url=json_url, headers=headers, auth=igvf_auth).json()

    out_dict: Dict = {'analysis_set': analysis_accession}
    # loop through fields to get
    files = json_results['files']
    for file in files:
        if file['content_type'] == 'sparse gene count matrix':
            out_dict['counts_matrix'] = file['accession']
        elif file['content_type'] == 'fragments':
            out_dict['fragments'] = file['accession']

    with outfile.open('w') as outf:
        json.dump(out_dict, outf, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make jsons with accession IDs for downloading "
                                                 "the counts matrix (scRNA-seq) and fragments file"
                                                 "(scATAC-seq) for a given analysis accession ID")
    parser.add_argument("--access_key", type=str, help="Access key from IGVF portal", required=True)
    parser.add_argument("--secret_key", type=str, help="Secret key from IGVF portal", required=True)
    parser.add_argument("-a", "--accession", type=str, help="Analysis accession")
    parser.add_argument("-o", "--outfile", type=Path, help="Where to save output")
    parser.add_argument("-u", "--base_url", type=str, help="base url", required=False,
                        default="https://api.data.igvf.org/analysis-sets")

    args = parser.parse_args()

    get_igvf_ids_for_download(
        analysis_accession=args.accession,
        access_key=args.access_key,
        secret_key=args.secret_key,
        base_url=args.base_url,
        outfile=args.outfile
    )

