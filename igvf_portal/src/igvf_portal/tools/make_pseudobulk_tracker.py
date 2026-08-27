import csv
import logging
import os
from collections.abc import Iterable, Iterator, Mapping
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import MappingProxyType
from typing import cast

import gspread
from igvf_client import IgvfApi
from igvf_client.models.analysis_set import AnalysisSet
from igvf_client.models.pseudobulk_set import PseudobulkSet
from igvf_client.models.search_result_item import SearchResultItem
from igvf_client.models.tabular_file import TabularFile

from igvf_portal import utils
from igvf_portal.annotations_file_qc import AnnotationsFileQc
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode, PseudobulkUploadStatus
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId, PseudobulkTrackerRow

_STATUS_ROW_COLORS: Mapping[PseudobulkUploadStatus, Mapping[str, float]] = (
    MappingProxyType(
        {
            PseudobulkUploadStatus.NEEDS_FIX: MappingProxyType(
                {"red": 0.957, "green": 0.8, "blue": 0.8}
            ),
            PseudobulkUploadStatus.UNATTEMPTED: MappingProxyType(
                {"red": 1.0, "green": 0.949, "blue": 0.8}
            ),
            PseudobulkUploadStatus.COMPLETE: MappingProxyType(
                {"red": 0.851, "green": 0.918, "blue": 0.827}
            ),
        }
    )
)


def _clear_conditional_format_rules(
    spreadsheet: gspread.Spreadsheet, worksheet: gspread.Worksheet
) -> None:
    """Remove any pre-existing conditional format rules on the worksheet so they don't accumulate on re-upload."""
    metadata = spreadsheet.fetch_sheet_metadata()
    sheet_properties = next(
        sheet
        for sheet in metadata["sheets"]
        if sheet["properties"]["sheetId"] == worksheet.id
    )
    num_rules = len(sheet_properties.get("conditionalFormats", []))
    if num_rules == 0:
        return
    # Delete from the highest index down, since deleting shifts later indices.
    requests = [
        {"deleteConditionalFormatRule": {"sheetId": worksheet.id, "index": index}}
        for index in reversed(range(num_rules))
    ]
    spreadsheet.batch_update({"requests": requests})


def _add_status_conditional_formatting(
    spreadsheet: gspread.Spreadsheet,
    worksheet: gspread.Worksheet,
    num_rows: int,
    num_cols: int,
    status_col_index: int,
) -> None:
    """Color each data row by its "pseudobulking status" value."""
    status_col_letter = gspread.utils.rowcol_to_a1(1, status_col_index + 1).rstrip("1")
    grid_range = {
        "sheetId": worksheet.id,
        "startRowIndex": 1,  # Skip the header row.
        "endRowIndex": num_rows,
        "startColumnIndex": 0,
        "endColumnIndex": num_cols,
    }
    requests = [
        {
            "addConditionalFormatRule": {
                "rule": {
                    "ranges": [grid_range],
                    "booleanRule": {
                        "condition": {
                            "type": "CUSTOM_FORMULA",
                            "values": [
                                {
                                    "userEnteredValue": f'=${status_col_letter}2="{status.value}"'
                                }
                            ],
                        },
                        "format": {"backgroundColor": {**color}},
                    },
                },
                "index": index,
            }
        }
        for index, (status, color) in enumerate(_STATUS_ROW_COLORS.items())
    ]
    spreadsheet.batch_update({"requests": requests})


def upload_dataframe_to_google_sheet(
    rows_tsv: Path,
    spreadsheet_url: str,
    worksheet_name: str,
    credentials_path: Path,
) -> None:
    """Upload the rows in the local TSV to the specified spreadsheet URL."""
    # Authenticate using the service-account JSON file.
    client = gspread.service_account(filename=credentials_path)

    # Open the existing Google Sheet.
    spreadsheet = client.open_by_url(spreadsheet_url)

    rows = utils.read_tsv(rows_tsv)
    fields = list(PseudobulkTrackerRow.__annotations__.keys())
    num_cols = len(fields)
    num_rows = len(rows)
    # Use the requested worksheet, creating it if necessary.
    try:
        worksheet = spreadsheet.worksheet(worksheet_name)
    except gspread.WorksheetNotFound:
        worksheet = spreadsheet.add_worksheet(
            title=worksheet_name,
            rows=max(num_rows, 100),
            cols=max(num_cols, 10),
        )

    # Replace the previous contents.
    worksheet.clear()
    worksheet.update(rows, f"A1:{gspread.utils.rowcol_to_a1(num_rows, num_cols)}")

    _clear_conditional_format_rules(spreadsheet, worksheet)
    _add_status_conditional_formatting(
        spreadsheet,
        worksheet,
        num_rows=num_rows,
        num_cols=num_cols,
        status_col_index=fields.index("pseudobulking status"),
    )


def _write_tsv(outfile: Path, rows: Iterable[PseudobulkTrackerRow]) -> None:
    """Write rows as a TSV at the specified path."""
    outfile.parent.mkdir(exist_ok=True)
    with outfile.open("wt") as f_out:
        writer = csv.DictWriter(
            f_out,
            fieldnames=PseudobulkTrackerRow.__annotations__.keys(),
            delimiter="\t",
            quoting=csv.QUOTE_NONE,
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            f_out.flush()


def _get_pseudobulk_status_and_date(
    api: IgvfApi, igvf_lookup: IgvfLookup, analysis_set: AnalysisSet
) -> tuple[PseudobulkUploadStatus, str | None]:
    """Given an input AnalysisSet record, identify status of any PseudobulkSets and date of upload."""
    if analysis_set.input_for is None:
        return PseudobulkUploadStatus.UNATTEMPTED, None

    audit_failures = set()
    date = None
    for input_id in analysis_set.input_for:
        input_record = utils.retry(num_tries=3)(api.get_by_id)(input_id).actual_instance
        if input_record is None:
            raise RuntimeError(f"Unable to get actual_instance of {input_id}")
        if input_record.type is not None and "PseudobulkSet" in input_record.type:
            input_record = cast(PseudobulkSet, input_record)
            if input_record.creation_timestamp is not None:
                date = input_record.creation_timestamp
                accession = input_record.accession
                if accession is None:
                    raise RuntimeError("Unable to get accession from record")
            audit_failures.update(
                igvf_lookup.lookup_record(AccessionId(accession))["audit"].keys()
            )
            if len(audit_failures) > 0:
                break
    match len(audit_failures), date:
        case _, None:
            return PseudobulkUploadStatus.UNATTEMPTED, None
        case 0, _:
            return PseudobulkUploadStatus.COMPLETE, date
        case _:
            return PseudobulkUploadStatus.NEEDS_FIX, date


def _process_analysis_set(
    api: IgvfApi, igvf_lookup: IgvfLookup, analysis_set: AnalysisSet
) -> tuple[TabularFile | None, PseudobulkUploadStatus, str | None]:
    """Get annotations file, pseudobulk upload status, and upload date for the provided AnalysisSet record."""
    annotations_file: TabularFile | None = None
    if analysis_set.files is None:
        return None, PseudobulkUploadStatus.UNATTEMPTED, None
    for file_id in analysis_set.files:
        file = utils.retry(num_tries=3)(api.get_by_id)(file_id)
        if file.content_type == "cell annotations":
            annotations_file = cast(TabularFile | None, file.actual_instance)
            if annotations_file is not None and annotations_file.status == "deleted":
                annotations_file = None
            break
    if annotations_file is None:
        return None, PseudobulkUploadStatus.UNATTEMPTED, None
    pseudobulk_status, upload_date = _get_pseudobulk_status_and_date(
        api, igvf_lookup, analysis_set
    )
    return annotations_file, pseudobulk_status, upload_date


def _get_tracker_row(
    search_result_item: SearchResultItem,
) -> tuple[str | None, PseudobulkTrackerRow | None, frozenset[str] | None]:
    # Get analysis set
    analysis_set = cast(AnalysisSet, search_result_item.actual_instance)
    annotations_file, pseudobulk_status, upload_date = _process_analysis_set(
        api, igvf_lookup, analysis_set
    )

    # Save annotation file info
    if annotations_file is None:
        return analysis_set.accession, None, None
    else:
        annotations_file_qc = AnnotationsFileQc.qc(api, annotations_file.accession)
        row: PseudobulkTrackerRow = {
            "principal analysis set accession": f"{analysis_set.accession}",
            "annotation file accession": f"{annotations_file.accession}",
            "lab": f"{annotations_file.lab}",
            "pseudobulking status": pseudobulk_status.value,
            "processed date": "" if upload_date is None else upload_date,
            "missing annotations columns": ",".join(
                annotations_file_qc.missing_columns
            ),
            "uniform pipeline status": annotations_file_qc.uniform_pipeline_status,
        }
        return analysis_set.accession, row, annotations_file_qc.cl_ids


def _set_num_workers(requested_workers: int | None, default: int = 12) -> int:
    if requested_workers is None or requested_workers <= 0:
        num_workers = os.process_cpu_count()
        if num_workers is None:
            num_workers = os.cpu_count()
            if num_workers is None:
                num_workers = default
        return num_workers
    else:
        return requested_workers


def _iter_rows(
    search_result_items: Iterable[SearchResultItem] | None,
    logger: logging.Logger,
    igvf_mode: IgvfMode,
    cl_ids_to_update: set[str],
    num_workers: int,
) -> Iterator[PseudobulkTrackerRow]:
    """Iterate over principal analysis sets and yield a PseudobulkTrackerRow with details of their status."""
    if search_result_items is None:
        raise ValueError("No principal analysis sets were found")
    logger.info(f"Scanning analysis sets with {num_workers} workers.")
    with ProcessPoolExecutor(
        max_workers=num_workers, initializer=_worker_initalizer, initargs=(igvf_mode,)
    ) as executor:
        for result in executor.map(
            _get_tracker_row, search_result_items, buffersize=2 * num_workers
        ):
            match result:
                case accession, None, _:
                    logger.warning(
                        f"Found analysis set {accession} with no annotation file."
                    )
                case accession, row, cl_ids:
                    date_info = (
                        ""
                        if len(row["processed date"]) == 0
                        else f" (uploaded {row['processed date']})"
                    )
                    num_cl_ids = 0 if cl_ids is None else len(cl_ids)
                    logger.info(
                        f"Found analysis set {accession} with status '{row['pseudobulking status']}'{date_info} and {num_cl_ids} unique CL_ids."
                    )
                    if cl_ids is not None:
                        cl_ids_to_update.update(cl_ids)
                    yield row


api: IgvfApi
igvf_lookup: IgvfLookup


def _worker_initalizer(igvf_mode: IgvfMode):
    global api
    global igvf_lookup
    utils.fix_igvf_logging(level=logging.WARNING)
    api = utils.open_igvf_api(igvf_mode=igvf_mode)
    igvf_lookup = IgvfLookup.new(igvf_mode=igvf_mode)


def make_pseudobulk_tracker(
    *,
    output: Path = Path("pseudobulk_status.tsv"),
    cli_ids_out: Path = Path("missing_cl_ids.txt"),
    compute: bool = True,
    upload: bool = False,
    spreadsheet_url: str = "https://docs.google.com/spreadsheets/d/1pYhl_ffq3pz_Y9zXmhqU62Ob5cLmzuQo_l5r6Bm7bZ0/edit",
    google_service_account_json: Path | None = None,
    num_workers: int = 12,
    igvf_mode: IgvfMode = IgvfMode.prod,
) -> None:
    """Write pseudobulk data set status to a table.

    To upload to google sheets, yout must pass --upload, as well as specify spreadsheet_url and google_service_account_json.

    Args:
        output: Path to write output table (as a TSV).
        compute: If True, compute statuses and write to TSV. If False, use existing local TSV for upload.
        upload: If True, upload to google sheets.
        spreedsheet_url: URL to google sheet to update.
        google_service_account_json: Path to google service account JSON.
    """
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")
    if upload and (
        google_service_account_json is None or not google_service_account_json.exists()
    ):
        raise ValueError(
            "google_service_account_json must point to a valid credentials JSON to upload."
        )

    if compute:
        # Log in to IGVF Portal
        api = utils.open_igvf_api(igvf_mode=igvf_mode)

        # Query all primary analysis sets with cell annotations
        search_results = api.search(
            type=["AnalysisSet"],
            limit="all",
            field_filters={
                "status!": "deleted",
                "file_set_type": "principal analysis",
                "files.content_type": "cell annotations",
            },
        )

        # Iterate through all queried primary analysis sets.
        # As a side effect update a set with all the CL_ids
        cl_ids = set()
        rows_iter = _iter_rows(
            search_results.graph,
            logger=logger,
            igvf_mode=igvf_mode,
            cl_ids_to_update=cl_ids,
            num_workers=_set_num_workers(requested_workers=num_workers, default=12),
        )
        _write_tsv(output, rows_iter)
        # now we've iterated through the rows and have all the unique CL_ids. Find the existing CL_ids on the portal
        search_results = api.search(
            type=["SampleTerm"], limit="all", field_filters={"status!": "deleted"}
        )
        portal_cl_ids = (
            set()
            if search_results.graph is None
            else {item.term_id for item in search_results.graph}
        )
        missing_ids = sorted(cl_ids.difference(portal_cl_ids))
        cli_ids_out.parent.mkdir(exist_ok=True, parents=True)
        with cli_ids_out.open("wt") as f_out:
            for missing_id in missing_ids:
                f_out.write(f"{missing_id}\n")

    if upload:
        if (
            google_service_account_json is None
            or not google_service_account_json.exists()
        ):
            raise ValueError(
                "google_service_account_json must point to a valid credentials JSON to upload."
            )
        upload_dataframe_to_google_sheet(
            output,
            spreadsheet_url=spreadsheet_url,
            worksheet_name="IGVF annotation files",
            credentials_path=google_service_account_json,
        )
