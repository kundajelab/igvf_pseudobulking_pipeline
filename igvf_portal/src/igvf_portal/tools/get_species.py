from collections.abc import Iterable, Iterator
from typing import Literal, cast

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import ContentType, IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId


def _get_reference_file_accessions(
    igvf_lookup: IgvfLookup, check_accessions: Iterable[AccessionId]
) -> Iterator[AccessionId]:
    check_types = frozenset({ContentType.FRAGMENTS.value, ContentType.MATRIX.value})
    all_check_accessions = set(check_accessions)
    remaining = all_check_accessions.copy()
    while len(remaining) > 0:
        check_accession = remaining.pop()
        check_record = igvf_lookup.lookup_record(check_accession)
        if "AnalysisSet" in check_record["@type"]:
            new_check = {
                file_record["accession"]
                for file_record in check_record["files"]
                if file_record["content_type"] in check_types
            }.difference(all_check_accessions)
            all_check_accessions.update(new_check)
            remaining.update(new_check)
        else:
            if check_record["content_type"] in check_types:
                yield from cast(list[AccessionId], check_record["reference_files"])


def _get_reference_file_species(
    igvf_lookup: IgvfLookup, reference_file_accessions: Iterable[AccessionId]
) -> set[Literal["human", "mouse"]]:
    reference_file_summaries = {
        igvf_lookup.lookup_record(reference_file_accession)["file_set"]["summary"]
        for reference_file_accession in reference_file_accessions
    }
    return {
        "human" if "sapiens" in summary else "mouse"
        for summary in reference_file_summaries
    }


def get_species(
    key: str,
    *,
    igvf_mode: IgvfMode = IgvfMode.prod,
) -> None:
    """From list of

    Form gene_info.csv and tss.tsv.


    Args:
        key: Comma-separated list of alias or accession ID to look for reference file dependencies
        igvf_mode: Mode for accessing the IGVF Portal.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)
    split_keys = {AccessionId(_split_key.strip()) for _split_key in key.split(",")}
    reference_file_accessions = set(
        _get_reference_file_accessions(
            igvf_lookup=igvf_lookup, check_accessions=split_keys
        )
    )
    species = _get_reference_file_species(
        igvf_lookup=igvf_lookup, reference_file_accessions=reference_file_accessions
    )
    match len(species):
        case 0:
            raise ValueError("Unable to find species description.")
        case 1:
            print(species.pop())
        case _:
            raise ValueError(f"Found multiple species: {', '.join(species)}")
