from collections.abc import Iterable, Iterator, Mapping
from types import MappingProxyType
from typing import Final, Literal

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import ContentType, IgvfMode
from igvf_portal.types import AccessionId, PortalId


def _get_reference_file_ids(
    connection: PConnection, check_accessions: Iterable[AccessionId]
) -> Iterator[PortalId]:
    check_types = frozenset({ContentType.FRAGMENTS, ContentType.MATRIX})
    all_check_accessions = set(check_accessions)
    remaining = all_check_accessions.copy()
    while len(remaining) > 0:
        check_accession = remaining.pop()
        check_record = connection.lookup_record(check_accession)
        if "AnalysisSet" in check_record["@type"]:
            new_check = {
                file_record["accession"]
                for file_record in check_record["files"]
                if any(content_type.is_wanted(file_record) for content_type in check_types)
            }.difference(all_check_accessions)
            all_check_accessions.update(new_check)
            remaining.update(new_check)
        else:
            if any(content_type.is_wanted(check_record) for content_type in check_types):
                for ref in check_record["reference_files"]:
                    yield ref if isinstance(ref, str) else ref["@id"]


_SPECIES_NAMES: Final[Mapping[Literal["human", "mouse"], str]] = MappingProxyType(
    {"human": "Homo sapiens", "mouse": "Mus musculus"}
)
"""The scientific name that a reference file set's summary gives for each supported species."""


def _get_summary_species(summary: str) -> Literal["human", "mouse"]:
    """Get the species named by a reference file set's summary, e.g. "Homo sapiens GRCh38 genome".

    Raises:
        ValueError: The summary names no supported species, or more than one. For example, the
            summary of a mixed human-mouse reference names neither: "GRCh38, mm10 ... genome".
    """
    named = [species for species, name in _SPECIES_NAMES.items() if name in summary]
    if len(named) != 1:
        raise ValueError(
            f"Unable to determine species from reference summary '{summary}': expected exactly "
            f"one of {', '.join(_SPECIES_NAMES.values())}."
        )
    return named[0]


def _get_reference_file_species(
    connection: PConnection, reference_file_ids: Iterable[PortalId]
) -> set[Literal["human", "mouse"]]:
    reference_file_summaries = {
        connection.lookup_record(reference_file_accession)["file_set"]["summary"]
        for reference_file_accession in reference_file_ids
    }
    return {_get_summary_species(summary) for summary in reference_file_summaries}


def get_species(
    key: str,
    *,
    igvf_mode: IgvfMode = IgvfMode.prod,
) -> None:
    """Infer the species of the reference files used by the given records. Display to stdout.

    Prints "human" or "mouse". Raises an error if no species, or more than one, is found.

    Args:
        key: Comma-separated list of alias or accession ID to look for reference file
            dependencies.
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    connection = PConnection.new(igvf_mode=igvf_mode)
    split_keys = {AccessionId(_split_key.strip()) for _split_key in key.split(",")}
    reference_file_ids = set(
        _get_reference_file_ids(connection=connection, check_accessions=split_keys)
    )
    species = _get_reference_file_species(
        connection=connection, reference_file_ids=reference_file_ids
    )
    match len(species):
        case 0:
            raise ValueError("Unable to find species description.")
        case 1:
            print(species.pop())
        case _:
            raise ValueError(f"Found multiple species: {', '.join(species)}")
