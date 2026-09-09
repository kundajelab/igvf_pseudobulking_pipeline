from logging import Logger
from pathlib import Path
from types import MappingProxyType
from typing import Final, TextIO, cast

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import ContentType, IgvfMode, MultipleRecordsAction
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId, IgvfRecord

HEADERS: Final[MappingProxyType[str, str]] = MappingProxyType(
    {"accept": "application/json"}
)
STDIN: Final[Path] = Path("-")


def _get_content_type_record(
    analysis_set_record: IgvfRecord,
    content_type: ContentType,
    multiple_records_action: MultipleRecordsAction,
    igvf_lookup: IgvfLookup,
    logger: Logger,
) -> IgvfRecord | None:
    content_type_records = [
        file_record
        for file_record in analysis_set_record["files"]
        if file_record["status"] != "deleted"
        and file_record["content_type"] == content_type.value
    ]
    for record in content_type_records:
        logger.info(
            f"Found {content_type.value} file with accession {record['accession']}."
        )
    match len(content_type_records), multiple_records_action:
        case 0 | 1, _:
            pass
        case _, MultipleRecordsAction.RAISE:
            pass
        case _, MultipleRecordsAction.KEEP_FILTERED:
            logger.info("Keeping filtered records.")
            # need to look up full record to get filtered status
            content_type_records = [
                record
                for record in content_type_records
                if igvf_lookup.lookup_record(record["accession"]).get("filtered", False)
            ]
        case _, MultipleRecordsAction.KEEP_UNFILTERED:
            logger.info("Keeping unfiltered records.")
            content_type_records = [
                record
                for record in content_type_records
                if not igvf_lookup.lookup_record(record["accession"]).get(
                    "filtered", False
                )
            ]
    for record in content_type_records:
        logger.info(
            f"Kept {content_type.value} file with accession {record['accession']}."
        )
    match len(content_type_records):
        case 0:
            return None
        case 1:
            return content_type_records[0]
        case _:
            raise ValueError(
                f"Found {len(content_type_records)} {content_type.value} files for {analysis_set_record['accession']}"
            )


def _write_download_entry(
    record: IgvfRecord | None,
    content_type: ContentType,
    igvf_mode: IgvfMode,
    igvf_portal_region: str,
    accession: AccessionId,
    f_out: TextIO,
    logger: Logger,
) -> None:
    """Use the specified record to write an aria2c download entry."""
    if record is None:
        logger.info(f"Got no {content_type.value} files.")
    else:
        url = utils.get_record_http_url(
            record,
            igvf_portal_region=igvf_portal_region,
            igvf_mode=igvf_mode,
        )
        logger.info(f"\tGot URL for {content_type.value}: {url}")
        f_out.write(f"{url}\n")
        f_out.write(f"  out={accession}.{content_type.extension}\n")


def get_url(
    accession: str,
    *,
    igvf_portal_region: str = "us-west-2",
    output: Path = STDIN,
    igvf_mode: IgvfMode = IgvfMode.prod,
    multiple_records_action: MultipleRecordsAction = MultipleRecordsAction.KEEP_UNFILTERED,
    accession_delimiter: str = ";",
):
    """Get download URLs for raw RNA h5ad and fragments bed.gz.

    Three forms are acceptible input accessions:
    1. Input accession can be for an AnalysisSet, in which case the IGVF Portal is queried to find
    matrix and fragment files.
    2. Alternatively it may be the accession for a matrix or fragments file. That is downloaded to
    a file name determined by its file_set accession.
    3. It can be an `accession_delimiter`-separated tuple specifying input_accession and
    output_accession. This works like cases 1 and 2 except that the output file-name is set to use
    the output accession as its base.

    Use the IGVF portal API to get the download URLs (which will be downloaded by aria2c).
    IGVF_API_KEY and IGVF_SECRET_KEY environment variables must be set to successfully
    log into the IGVF portal API.

    Args:
        accession: Analysis set accession to get URLs for.
        igvf_portal_region: S3 region for the IGVF portal S3 backing.
        output: File to write URLs to. If "-", write to stdout.
        igvf_mode: "prod", "staging", or "sandbox"
        multiple_records_action: If more than one matrix/fragments files are found:
          if KEEP_FILTERED, keep filtered records
          if KEEP_UNFILTERED, keep unfiltered records
          if RAISE, raise an exception if there are multiple records
          If only one record is found, keep it regardless of this setting.
        accession_delimiter: use this character to try to split accessions into input/output tuples
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    igvf_lookup = IgvfLookup.new(igvf_mode)
    if accession_delimiter in accession:
        input_accession, output_accession = cast(
            tuple[AccessionId, AccessionId], accession.split(accession_delimiter, 1)
        )
    else:
        input_accession = AccessionId(accession)
        output_accession = None
    input_record = igvf_lookup.lookup_record(input_accession)

    if f"{output}" == "-":
        output = Path("/dev/stdout")

    with output.open("at") as f_out:
        if "AnalysisSet" in input_record["@type"]:
            if output_accession is None:
                output_accession = AccessionId(input_record["accession"])
            # lookup matrices and fragments from portal
            for content_type in (ContentType.MATRIX, ContentType.FRAGMENTS):
                record = _get_content_type_record(
                    analysis_set_record=input_record,
                    content_type=content_type,
                    multiple_records_action=multiple_records_action,
                    igvf_lookup=igvf_lookup,
                    logger=logger,
                )
                _write_download_entry(
                    record=record,
                    content_type=content_type,
                    igvf_mode=igvf_mode,
                    igvf_portal_region=igvf_portal_region,
                    accession=output_accession,
                    f_out=f_out,
                    logger=logger,
                )
        elif "content_type" in input_record:
            if input_record["content_type"] in (
                ContentType.MATRIX.value,
                ContentType.FRAGMENTS.value,
            ):
                if output_accession is None:
                    output_accession = input_record["file_set"]["accession"]
                _write_download_entry(
                    record=input_record,
                    content_type=ContentType[input_record["content_type"]],
                    igvf_mode=igvf_mode,
                    igvf_portal_region=igvf_portal_region,
                    accession=output_accession,
                    f_out=f_out,
                    logger=logger,
                )
            else:
                raise ValueError(
                    f"Got accession for record with content_type: {input_record['content_type']}"
                )
        else:
            raise ValueError(
                f"Got accession for record with @type: {input_record.get('@type')} that does not have a content_type."
            )
