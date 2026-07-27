from pathlib import Path

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import Alias


def download_file(
    key: str,
    *,
    igvf_mode: IgvfMode = IgvfMode.prod,
    output: Path | None = None,
    chunk_size: int = 8192,
) -> None:
    """Infer the principal analysis accession IDs from the input metadata file. Display to stdout

    Args:
        metadata_file: Path to annotations file.
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)
    record = igvf_lookup.lookup_record(Alias(key))
    utils.download_record(
        record=record,
        igvf_mode=igvf_mode,
        chunk_size=chunk_size,
        output=output,
        logger=logger,
    )
