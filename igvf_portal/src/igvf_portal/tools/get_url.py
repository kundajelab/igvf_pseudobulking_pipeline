from pathlib import Path
from types import MappingProxyType
from typing import Final

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import ContentType, IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId, IgvfRecord


def _http_url(file_record: IgvfRecord, igvf_portal_region: str) -> str:
    s3_uri = file_record["s3_uri"]
    _, _, bucket_name, object_key = s3_uri.split("/", 3)
    return f"https://{bucket_name}.s3.{igvf_portal_region}.amazonaws.com/{object_key}"


HEADERS: Final[MappingProxyType[str, str]] = MappingProxyType(
    {"accept": "application/json"}
)
STDIN: Final[Path] = Path("-")


def get_url(
    accession: str,
    *,
    igvf_portal_region: str = "us-west-2",
    output: Path = STDIN,
    igvf_mode: IgvfMode = IgvfMode.prod,
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
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    igvf_lookup = IgvfLookup.new(igvf_mode)
    analysis_set_record = igvf_lookup.lookup_record(AccessionId(accession))

    if f"{output}" == "-":
        output = Path("/dev/stdout")

    with output.open("at") as f_out:
        num_matrices = 0
        num_fragments = 0
        for file_record in analysis_set_record["files"]:
            if file_record["status"] == "deleted":
                continue
            match file_record["content_type"]:
                case ContentType.MATRIX.value:
                    num_matrices += 1
                    logger.info(f"Got {ContentType.MATRIX.value} for {accession}")
                    url = _http_url(file_record, igvf_portal_region=igvf_portal_region)
                    logger.info(f"\tGot URL for {ContentType.MATRIX.value}: {url}")
                    f_out.write(f"{url}\n")
                    f_out.write(f"  out={accession}.h5ad\n")
                case ContentType.FRAGMENTS.value:
                    num_fragments += 1
                    logger.info(f"Got {ContentType.FRAGMENTS.value} for {accession}")
                    url = _http_url(file_record, igvf_portal_region=igvf_portal_region)
                    logger.info(f"\tGot URL for {ContentType.MATRIX.value}: {url}")
                    f_out.write(f"{url}\n")
                    f_out.write(f"  out={accession}.bed.gz\n")
                case _:
                    logger.info(
                        f"Got unused content type for {accession}: "
                        f"{file_record['content_type']}"
                    )
        if num_fragments == num_matrices == 0:
            raise ValueError(f"No data found for {accession}")
        if num_fragments > 1:
            raise ValueError(
                f"Found {num_fragments} {ContentType.FRAGMENTS.value} files for {accession}"
            )
        if num_matrices > 1:
            raise ValueError(
                f"Found {num_matrices} {ContentType.MATRIX.value} files for {accession}"
            )
