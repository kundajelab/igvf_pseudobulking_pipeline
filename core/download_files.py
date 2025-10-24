#!/usr/bin/python3
# Copyright 2025 Diane Trout - Caltech
# SPDX-License-Identifier: BSD-3-Clause

from argparse import ArgumentParser
import logging
from netrc import netrc
import os
from pathlib import Path
import requests
import shutil
from urllib.parse import urlparse, urljoin


LOGGER = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)

    auth = get_igvf_auth(args.base_url, args.access_key, args.secret_key)
    outdir = Path(args.outdir)

    for accession in args.accessions:
        file_info = get_file_info(args.base_url, accession, auth)
        if "counts_matrix_href" in file_info:
            download_rna_counts_matrix(outdir, file_info, auth, args.dry_run)
        if "fragments_href" in file_info:
            download_atac_fragments(outdir, file_info, auth, args.dry_run)


def make_parser():
    parser = ArgumentParser()
    parser.add_argument("accessions", nargs="+", help="accessions to scan")
    parser.add_argument("-o", "--outdir", help="target directory", required=True)
    parser.add_argument(
        "-u",
        "--base_url",
        type=str,
        help="base url",
        required=False,
        default="https://api.data.igvf.org/analysis-sets",
    )
    parser.add_argument("--access_key", type=str, help="Access key from IGVF portal")
    parser.add_argument("--secret_key", type=str, help="Secret key from IGVF portal")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would be done without downloading",
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="verbose level logging"
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="debug level logging"
    )
    return parser


def get_igvf_auth(base_url, user=None, secret=None):
    """Get the authentication tokens from several possible locations

    It checks the arguments passed in, the IGVF_API_KEY/IGVF_SECRET_KEY
    environment variables, and the netrc file

    A netrc file has the format
    machine HOSTNAME login ACCESS_KEY password SECRET_KEY
    """
    if user is not None and secret is not None:
        return (user, secret)

    user = os.environ.get("IGVF_API_KEY", None)
    secret = os.environ.get("IGVF_SECRET_KEY", None)

    if user is not None and secret is not None:
        return (user, secret)

    # it might be fancy to include keyring sources too in here

    parts = urlparse(base_url)
    nc = netrc()
    auth_tokens = nc.authenticators(parts.netloc)
    if auth_tokens is not None:
        return (auth_tokens[0], auth_tokens[2])

    raise ValueError("Need to specify access and secret key")


def get_file_info(base_url, accession, auth):
    fragment = f"analysis-sets/{accession}/"
    analysis_set_url = urljoin(base_url, fragment)
    LOGGER.info("Downloading analysis set {} metadata".format(analysis_set_url))
    resp = requests.get(analysis_set_url, auth=auth)

    if resp.status_code != 200:
        LOGGER.error("HTTP Error: {} {}".format(resp.status_code, resp.text))
        return None

    data = resp.json()
    results = {
        "analysis_accession": accession,
        "assay_titles": data.get("assay_titles"),
    }

    for f in data.get("files", []):
        content_type = f["content_type"]
        if content_type == "sparse gene count matrix":
            results["counts_matrix_accession"] = f["accession"]
            results["counts_matrix_href"] = urljoin(base_url, f["href"])
        elif content_type == "fragments":
            results["fragments_accession"] = f["accession"]
            results["fragments_href"] = urljoin(base_url, f["href"])

    return results


def download_rna_counts_matrix(outdir, file_info, auth, dry_run):
    target_dir = Path(outdir) / "raw_rna"
    if not target_dir.exists():
        target_dir.mkdir(parents=True)

    target_name = target_dir / (file_info["analysis_accession"] + ".h5ad")
    download_file(file_info["counts_matrix_href"], target_name, auth, dry_run)


def download_atac_fragments(outdir, file_info, auth, dry_run):
    target_dir = Path(outdir) / "raw_fragments"
    if not target_dir.exists():
        target_dir.mkdir(parents=True)

    target_name = target_dir / (file_info["analysis_accession"] + ".bed.gz")
    download_file(file_info["fragments_href"], target_name, auth, dry_run)


def download_file(href, target_filename, auth, dry_run):
    LOGGER.debug("Downloading {} to {}".format(href, target_filename))

    if not dry_run:
        resp = requests.get(href, stream=True, auth=auth)
        if resp.status_code == 200:
            with open(target_filename, "wb") as outstream:
                shutil.copyfileobj(resp.raw, outstream)
        else:
            LOGGER.error("HTTP Error: {} {}".format(resp.status_code, resp.text))


if __name__ == "__main__":
    main()
