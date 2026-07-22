from pathlib import Path
from types import MappingProxyType
from typing import Final, cast
import logging
import os
import requests

import defopt

HEADERS: Final[MappingProxyType[str, str]] = MappingProxyType({"accept": "application/json"})



def get_url(
    accession: str,
    *,
    igvf_portal_url: str = "https://api.data.igvf.org",
    analysis_sets_folder: str = "analysis-sets",
    output: Path = Path("-"),
    get_redirect: bool = True,
):
    """Get download URLs for raw RNA h5ad and fragments bed.gz for a given analysis set accession.

    Use the IGVF portal API to get the download URLs (which will be downloaded by aria2c).
    IGVF_API_KEY and IGVF_SECRET_KEY environment variables must be set to successfully
    log into the IGVF portal API.

    Args:
        accession: Analysis set accession to get URLs for.
        igvf_portal_url: Base URL to query IGVF portal.
        analysis_sets_folder: Subfolder listing analysis sets.
        output: File to write URLs to. If "-", write to stdout.
        get_redirect: If True, get the redirect URL at S3 (public but with an expiration time),
            if False, just get the IGVF portal download URL
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger("get-url")

    igvf_auth: tuple[str, str] = (
        os.environ["IGVF_API_KEY"],
        os.environ["IGVF_SECRET_KEY"],
    )
    json_url = f"{igvf_portal_url}/{analysis_sets_folder}/{accession}/?format=json"
    json_results = requests.get(url=json_url, headers=HEADERS, auth=igvf_auth).json()

    if f"{output}" == "-":
        output = Path("/dev/stdout")

    with output.open("at") as f_out:
        for json_result_file in cast(list[dict[str, str]], json_results["files"]):
            match json_result_file["content_type"]:
                case "sparse gene count matrix" | "cell by gene matrix":
                    logger.info(f"Got sparse gene count matrix for {accession}")
                    if get_redirect:
                        counts_matrix_url: str = json_result_file["s3_uri"]
                        logger.info(f"Got redirect URL for counts matrix: {counts_matrix_url}")
                    else:
                        counts_matrix_url = f"{igvf_portal_url}{json_result_file['href']}"
                        logger.info(f"Got IGVF portal URL: {igvf_portal_url}")
                    f_out.write(f"{counts_matrix_url}\n")
                    f_out.write(f"  out={accession}.h5ad\n")
                case "fragments":
                    logger.info(f"Got fragments for {accession}")
                    if get_redirect:
                        fragments_url: str = json_result_file["s3_uri"]
                        logger.info(f"Got redirect URL for fragments: {fragments_url}")
                    else:
                        fragments_url = f"{igvf_portal_url}{json_result_file['href']}"
                        logger.info(f"Got IGVF portal URL: {igvf_portal_url}")
                    f_out.write(f"{fragments_url}\n")
                    f_out.write(f"  out={accession}.bed.gz\n")
                case _:
                    logger.info(
                        f"Got unused content type for {accession}: "
                        f"{json_result_file['content_type']}"
                    )


def main():
    defopt.run(get_url)
