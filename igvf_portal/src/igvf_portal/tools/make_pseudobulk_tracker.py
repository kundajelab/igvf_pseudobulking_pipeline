import csv
import gzip
import logging
import os
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, wait
from io import StringIO
from pathlib import Path
from typing import TextIO, cast

import gspread
import igvf_utils
from igvf_client import ApiClient, Configuration, IgvfApi
from igvf_client.models.analysis_set import AnalysisSet
from igvf_client.models.pseudobulk_set import PseudobulkSet
from igvf_client.models.search_result_item import SearchResultItem
from igvf_client.models.tabular_file import TabularFile

from igvf_portal import utils
from igvf_portal.constants import VERSION
from igvf_portal.enums import PseudobulkUploadStatus
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import AccessionId, PseudobulkTrackerRow

_STATUS_ROW_COLORS: dict[PseudobulkUploadStatus, dict[str, float]] = {
    PseudobulkUploadStatus.NEEDS_FIX: {"red": 0.957, "green": 0.8, "blue": 0.8},
    PseudobulkUploadStatus.UNATTEMPTED: {"red": 1.0, "green": 0.949, "blue": 0.8},
    PseudobulkUploadStatus.COMPLETE: {"red": 0.851, "green": 0.918, "blue": 0.827},
}

_REQUIRED_ANNOTATION_COLUMNS: frozenset[str] = frozenset(
    {
        "barcode_sample",
        "cell_name",
        "cell_description",
        "CL_id",
        "CL_term_name",
        "subsample",
        "analysis_set_accession",
    }
)


def _read_tsv_bytes(f_in: TextIO) -> list[list[str]]:
    reader = csv.reader(f_in, delimiter="\t")
    return [row for row in reader]


def _read_tsv(tsv: Path) -> list[list[str]]:
    with tsv.open("rt") as f_in:
        return _read_tsv_bytes(f_in)


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
                        "format": {"backgroundColor": color},
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

    rows = _read_tsv(rows_tsv)
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
        input_record = api.get_by_id(input_id).actual_instance
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


def _collect_uniform_pipeline_status(
    analysis_sets: Iterable[AnalysisSet],
) -> set[str] | None:
    statuses = set()
    for analysis_set in analysis_sets:
        status = analysis_set.uniform_pipeline_status
        match status:
            case None:
                return None
            case _:
                statuses.add(status)
    return statuses


def _qc_annotations_file(
    api: IgvfApi, annotations_file: TabularFile
) -> tuple[set[str], set[str] | None]:
    if annotations_file.accession is None:
        raise ValueError("Missing accession for annotations file")
    annotations_compressed_bytes = api.download(annotations_file.accession)
    annotations_decompressed_str = gzip.decompress(
        annotations_compressed_bytes
    ).decode()
    tsv_rows = _read_tsv_bytes(StringIO(annotations_decompressed_str))
    missing_columns = set(_REQUIRED_ANNOTATION_COLUMNS.difference(tsv_rows[0]))

    if "analysis_set_accession" in missing_columns:
        on_uniformly_processed_data = None
    else:
        column_index = next(
            idx
            for idx, column in enumerate(tsv_rows[0])
            if column == "analysis_set_accession"
        )
        analysis_set_accessions = {row[column_index] for row in tsv_rows[1:]}
        on_uniformly_processed_data = _collect_uniform_pipeline_status(
            cast(AnalysisSet, api.get_by_id(x).actual_instance)
            for x in analysis_set_accessions
        )
    return missing_columns, on_uniformly_processed_data


def _process_analysis_set(
    api: IgvfApi, igvf_lookup: IgvfLookup, analysis_set: AnalysisSet
) -> tuple[TabularFile | None, PseudobulkUploadStatus, str | None]:
    """Get annotations file, pseudobulk upload status, and upload date for the provided AnalysisSet record."""
    annotations_file: TabularFile | None = None
    if analysis_set.files is None:
        return None, PseudobulkUploadStatus.UNATTEMPTED, None
    for file_id in analysis_set.files:
        file = api.get_by_id(file_id)
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
) -> tuple[str | None, PseudobulkTrackerRow | None]:
    # Get analysis set
    analysis_set = cast(AnalysisSet, search_result_item.actual_instance)
    annotations_file, pseudobulk_status, upload_date = _process_analysis_set(
        api, igvf_lookup, analysis_set
    )

    # Save annotation file info
    if annotations_file is None:
        return analysis_set.accession, None
    else:
        missing_columns, uniform_pipeline_status = _qc_annotations_file(
            api, annotations_file
        )
        row: PseudobulkTrackerRow = {
            "principal analysis set accession": f"{analysis_set.accession}",
            "annotation file accession": f"{annotations_file.accession}",
            "lab": f"{annotations_file.lab}",
            "pseudobulking status": pseudobulk_status.value,
            "processed date": "" if upload_date is None else upload_date,
            "missing annotations columns": ",".join(missing_columns),
            "uniform pipeline status": "UNKNOWN"
            if uniform_pipeline_status is None
            else ",".join(uniform_pipeline_status),
        }
        return analysis_set.accession, row


def _iter_rows(
    search_result_items: Iterable[SearchResultItem] | None,
    logger: logging.Logger,
    num_workers: int | None = None,
) -> Iterator[PseudobulkTrackerRow]:
    """Iterate over principal analysis sets and yield a PseudobulkTrackerRow with details of their status."""
    if search_result_items is None:
        raise ValueError("No principal analysis sets were found")
    if num_workers is None or num_workers <= 0:
        num_workers = os.process_cpu_count()
        if num_workers is None:
            num_workers = os.cpu_count()
    with ProcessPoolExecutor(
        max_workers=num_workers, initializer=_worker_initalizer
    ) as executor:
        futures = [
            executor.submit(_get_tracker_row, item) for item in search_result_items
        ]
        for future in futures:
            wait([future])
            match future.result():
                case accession, None:
                    logger.warning(
                        f"Found analysis set {accession} with no annotation file."
                    )
                case accession, row:
                    date_info = (
                        ""
                        if len(row["processed date"]) == 0
                        else f" (uploaded {row['processed date']})"
                    )
                    logger.info(
                        f"Found analysis set {accession} with status '{row['pseudobulking status']}'{date_info}."
                    )
                    yield row

    # for search_result_item in search_result_items:
    #     # Get analysis set
    #     analysis_set = cast(AnalysisSet, search_result_item.actual_instance)
    #     annotations_file, pseudobulk_status, upload_date = _process_analysis_set(
    #         api, igvf_lookup, analysis_set
    #     )

    #     # Save annotation file info
    #     if annotations_file is None:
    #         logger.warning(
    #             f"Found analysis set {analysis_set.accession} with no annotation file."
    #         )
    #     else:
    #         missing_columns, uniform_pipeline_status = _qc_annotations_file(api, annotations_file)
    #         row: PseudobulkTrackerRow = {
    #             "principal analysis set accession": f"{analysis_set.accession}",
    #             "annotation file accession": f"{annotations_file.accession}",
    #             "lab": f"{annotations_file.lab}",
    #             "pseudobulking status": pseudobulk_status.value,
    #             "processed date": "" if upload_date is None else upload_date,
    #             "missing annotations columns": ",".join(missing_columns),
    #             "uniform pipeline status": "UNKNOWN" if uniform_pipeline_status is None else ",".join(uniform_pipeline_status)
    #         }
    #         date_info = "" if upload_date is None else f" (uploaded {upload_date})"
    #         logger.info(
    #             f"Found analysis set {analysis_set.accession} with status '{row['pseudobulking status']}'{date_info}."
    #         )
    #         yield row


def _open_api() -> tuple[IgvfApi, IgvfLookup]:
    utils.check_access_keys()
    config = Configuration(
        access_key=os.environ["IGVF_API_KEY"],
        secret_access_key=os.environ["IGVF_SECRET_KEY"],
    )
    client = ApiClient(config)
    api = IgvfApi(client)
    igvf_lookup = IgvfLookup.new("prod")
    return api, igvf_lookup


api: IgvfApi
igvf_lookup: IgvfLookup


def _worker_initalizer():
    global api
    global igvf_lookup
    utils.setup_logger(igvf_utils.debug_logger, logging.WARNING)
    utils.setup_logger(logging.getLogger(name="root"), logging.WARNING)
    api, igvf_lookup = _open_api()


def make_pseudobulk_tracker(
    *,
    output: Path = Path("pseudobulk_status.tsv"),
    compute: bool = True,
    upload: bool = False,
    spreadsheet_url: str = "https://docs.google.com/spreadsheets/d/1pYhl_ffq3pz_Y9zXmhqU62Ob5cLmzuQo_l5r6Bm7bZ0/edit",
    google_service_account_json: Path | None = None,
    num_workers: int = -1,
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
        api, _ = _open_api()

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

        # Iterate through all queried primary analysis sets
        rows_iter = _iter_rows(search_results.graph, logger, num_workers=num_workers)
        _write_tsv(output, rows_iter)

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
