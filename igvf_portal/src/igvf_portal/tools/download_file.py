from pathlib import Path

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode
from igvf_portal.types import Alias


def download_file(
    key: str,
    *,
    igvf_mode: IgvfMode = IgvfMode.prod,
    output: Path | None = None,
    chunk_size: int = 2**20,
) -> None:
    """Download the file described by an IGVF Portal record.

    Args:
        key: Alias or accession ID of the file record to download.
        igvf_mode: Mode for accessing the IGVF Portal.
        output: If specified, download to this path. If it is an existing folder or has no
            suffix, download into that folder using the record's href as file name. If
            unspecified, download in working folder.
        chunk_size: Chunk size for streaming download, in bytes.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    connection = PConnection.new(igvf_mode=igvf_mode)
    # downloading needs only the file's own href and s3_uri
    record = connection.lookup_record(Alias(key), frame="object")
    connection.download_record(record=record, chunk_size=chunk_size, output=output)
