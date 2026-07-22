from pathlib import Path

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId


def _get_accession_ids(metadata_file: Path) -> set[AccessionId]:
    """Read the accession IDs from the metadata file."""
    return {
        AccessionId(line["analysis_set_accession"])
        for line in utils.iter_csv_rows(metadata_file)
    }


def infer_principal_analysis(
    metadata_file: Path, igvf_mode: IgvfMode = IgvfMode.prod
) -> None:
    """Infer the principal analysis accession IDs from the input metadata file. Display to stdout

    Args:
        metadata_file: Path to annotations file.
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    accession_ids = _get_accession_ids(metadata_file)
    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)
    principal_accessions = igvf_lookup.infer_principal_accessions(accession_ids)
    print(",".join(sorted(principal_accessions)))
