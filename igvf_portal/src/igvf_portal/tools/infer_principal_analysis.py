import os
from pathlib import Path

from igvf_portal.utils import iter_csv_rows
from igvf_portal.enums import IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId


def _get_accession_ids(metadata_file: Path) -> set[AccessionId]:
    """Read the accession IDs from the metadata file."""
    return {
        AccessionId(line["analysis_set_accession"])
        for line in iter_csv_rows(metadata_file)
    }


def infer_principal_analysis(
    metadata_file: Path, igvf_mode: IgvfMode = IgvfMode.prod
) -> None:
    """Infer the principal analysis accession IDs from the input metadata file. Display to stdout

    Args:
        metadata_file: Path to annotations file.
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    if "IGVF_API_KEY" not in os.environ or "IGVF_SECRET_KEY" not in os.environ:
        raise ValueError("IGVF_API_KEY and IGVF_SECRET_KEY must be set in environment")
    accession_ids = _get_accession_ids(metadata_file)
    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)
    principal_accessions = igvf_lookup.infer_principal_accessions(accession_ids)
    print(",".join(sorted(principal_accessions)))
