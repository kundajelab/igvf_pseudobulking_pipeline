from pathlib import Path
from types import MappingProxyType
from typing import Final
import logging
import os
import requests

import defopt

HEADERS: Final[MappingProxyType[str, str]] = MappingProxyType({"accept": "application/json"})


def _get_igvf_portal_url(base_url: str, file_accession: str, extension: str) -> str:
    """Get the IGVF portal download URL for the specified data."""
    return f"{base_url}/{file_accession}/@@download/{file_accession}.{extension}"


def _get_redirect_url(igvf_portal_url: str, igvf_auth: tuple[str, str]) -> str:
    """Get the redirected public S3 URL for the desired file.

    Args:
        base_url: Base URL at IGVF portal for the desired file type.
        file_accession: Accession for this specific file (not analysis set accession)
        extension: output file extension
        igvf_auth: access key / secret key tuple to authorize access
    """
    redirect_response = requests.get(
        url=igvf_portal_url, headers=HEADERS, auth=igvf_auth, allow_redirects=False
    )
    if redirect_response.is_redirect:
        return redirect_response.headers.get("Location")
    else:
        raise ValueError(
            f"Expected a redirect response for URL {igvf_portal_url}, but got status code "
            f"{redirect_response.status_code}"
        )


def get_url(
    accession: str,
    *,
    igvf_metadata_url: str = "https://api.data.igvf.org/analysis-sets",
    igvf_fragments_url: str = "https://api.data.igvf.org/tabular-files",
    igvf_counts_matrix_url: str = "https://api.data.igvf.org/matrix-files",
    output: Path = Path("-"),
    get_redirect: bool = True,
):
    """Get download URLs for raw RNA h5ad and fragments bed.gz for a given analysis set accession.

    Use the IGVF portal API to get the download URLs (which will be downloaded by aria2c).
    IGVF_API_KEY and IGVF_SECRET_KEY environment variables must be set to successfully
    log into the IGVF portal API.

    Args:
        accession: Analysis set accession to get URLs for.
        igvf_metadata_url: Base URL to query IGVF portal for details of accession files.
        igvf_fragments_url: Base URL for IGVF fragments files.
        igvf_counts_matrix_url: Base URL for IGVF counts matrix files.
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
    json_url = f"{igvf_metadata_url}/{accession}/?format=json"
    json_results = requests.get(url=json_url, headers=HEADERS, auth=igvf_auth).json()

    if f"{output}" == "-":
        output = Path("/dev/stdout")

    with output.open("at") as f_out:
        for json_result_file in json_results["files"]:
            match json_result_file["content_type"]:
                case "sparse gene count matrix" | "cell by gene matrix":
                    logger.info(f"Got sparse gene count matrix for {accession}")
                    igvf_portal_url = _get_igvf_portal_url(
                        base_url=igvf_counts_matrix_url,
                        file_accession=json_result_file["accession"],
                        extension="h5ad",
                    )
                    if get_redirect:
                        counts_matrix_url = _get_redirect_url(igvf_portal_url, igvf_auth=igvf_auth)
                        logger.info(f"Got redirect URL for counts matrix: {counts_matrix_url}")
                    else:
                        counts_matrix_url = igvf_portal_url
                        logger.info(f"Got IGVF portal URL: {igvf_portal_url}")
                    f_out.write(f"{counts_matrix_url}\n")
                    f_out.write(f"  out={accession}.h5ad\n")
                case "fragments":
                    logger.info(f"Got fragments for {accession}")
                    igvf_portal_url = _get_igvf_portal_url(
                        base_url=igvf_fragments_url,
                        file_accession=json_result_file["accession"],
                        extension="bed.gz",
                    )
                    if get_redirect:
                        fragments_url = _get_redirect_url(igvf_portal_url, igvf_auth=igvf_auth)
                        logger.info(f"Got redirect URL for fragments: {fragments_url}")
                    else:
                        fragments_url = igvf_portal_url
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
