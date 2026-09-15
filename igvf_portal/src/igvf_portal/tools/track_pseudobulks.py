import csv
import logging
import os
from collections import defaultdict
from collections.abc import Iterable, Iterator, Mapping
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from threading import Lock, local
from types import MappingProxyType
from typing import Final, cast

import gspread
from igvf_client import IgvfApi
from igvf_client.models.analysis_set import AnalysisSet
from igvf_client.models.search_result_item import SearchResultItem
from igvf_client.models.tabular_file import TabularFile

from igvf_portal import utils
from igvf_portal.annotations_file_qc import AnnotationsFileQc
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import IgvfMode, PseudobulkUploadStatus
from igvf_portal.types import AccessionId, Alias, PortalId, PseudobulkTrackerRow

_THREAD_LOCAL_DATA: local = local()

_DEFAULT_NUM_WORKERS: Final[int] = 16

_STATUS_ROW_COLORS: Mapping[str, Mapping[str, float]] = MappingProxyType(
    {
        "red": MappingProxyType({"red": 0.99, "green": 0.6, "blue": 0.53}),
        "yellow": MappingProxyType({"red": 1.0, "green": 1.0, "blue": 0.63}),
        "green": MappingProxyType({"red": 0.18, "green": 0.85, "blue": 0.33}),
        "orange": MappingProxyType({"red": 0.99, "green": 0.85, "blue": 0.5333}),
    }
)

_PREFERRED_ASSAY_TITLES: frozenset[str] = frozenset(
    {
        "10x multiome with MULTI-seq",
        "10x multiome",
        "Parse SPLiT-seq",
        "snATAC-seq",
        "SHARE-seq",
        "10x snATAC-seq with Scale pre-indexing",
    }
)

_NO_PSEUDOBULKS: Final[str] = "NO_PSEUDOBULKS"


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
    uniform_pipeline_status_col_index: int,
) -> None:
    """Color data rows by both statuses, with failure conditions taking priority."""
    status_col_letter = gspread.utils.rowcol_to_a1(1, status_col_index + 1).rstrip("1")
    uniform_status_col_letter = gspread.utils.rowcol_to_a1(
        1, uniform_pipeline_status_col_index + 1
    ).rstrip("1")
    # Complete statuses have the form "complete: <num_pseudobulks>".
    is_complete = (
        f"REGEXMATCH(${status_col_letter}2,"
        rf'"^{PseudobulkUploadStatus.COMPLETE.value}: \d+$")'
    )
    # Sheets uses the first matching rule, so complete and not-validated statuses
    # take priority over red, including NO_ASSEMBLIES.
    rules = (
        (f"={is_complete}", _STATUS_ROW_COLORS["green"]),
        (
            f'=${status_col_letter}2="upload status not validated"',
            _STATUS_ROW_COLORS["orange"],
        ),
        (
            f"=OR(AND(${status_col_letter}2"
            f'<>"{PseudobulkUploadStatus.UNATTEMPTED.value}",NOT({is_complete})),'
            rf'REGEXMATCH(${uniform_status_col_letter}2,"\bNO_ASSEMBLIES\b"))',
            _STATUS_ROW_COLORS["red"],
        ),
        ("=TRUE", _STATUS_ROW_COLORS["yellow"]),
    )
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
                            "values": [{"userEnteredValue": formula}],
                        },
                        "format": {"backgroundColor": {**color}},
                    },
                },
                "index": index,
            }
        }
        for index, (formula, color) in enumerate(rules)
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
        uniform_pipeline_status_col_index=fields.index("uniform pipeline status"),
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


def _get_all_pseudobulk_audit_errors_and_timestamps(
    connection: PConnection,
    analysis_set: AnalysisSet,
    metadata_rows: list[list[str]],
) -> tuple[set[str], list[str]]:
    """Given an input AnalysisSet record, identify status of any PseudobulkSets and date of upload."""
    if analysis_set.input_for is None:
        return set(), []

    audit_failures: set[str] = set()
    submission_timestamps: list[str] = []

    for input_id in _get_input_ids(
        connection=connection, analysis_set=analysis_set, metadata_rows=metadata_rows
    ):
        _update_pseudobulk_audit_errors_and_date(
            input_id=PortalId(input_id),
            connection=connection,
            audit_failures=audit_failures,
            submission_timestamps=submission_timestamps,
        )

    return audit_failures, submission_timestamps


def _get_input_ids(
    connection: PConnection, analysis_set: AnalysisSet, metadata_rows: list[list[str]]
) -> set[PortalId]:
    to_check: set[PortalId] = (
        set()
        if analysis_set.input_for is None
        else cast(set[PortalId], set(analysis_set.input_for))
    )
    for intermediate_accession in AnnotationsFileQc.get_column_accessions(
        metadata_rows, "analysis_set_accession"
    ):
        to_check.update(connection.lookup_record(intermediate_accession)["input_for"])
    return to_check


def _update_pseudobulk_audit_errors_and_date(
    input_id: PortalId | AccessionId | Alias,
    connection: PConnection,
    audit_failures: set[str],
    submission_timestamps: list[str],
) -> None:
    input_record = connection.lookup_record(input_id)
    if "PseudobulkSet" not in input_record["@type"]:
        return

    submission_timestamps.append(input_record["submitted_files_timestamp"])
    audit_failures.update(
        f"{cast(dict[str, str], error).get('category')}"
        for error in cast(
            dict[str, dict[str, list[object]]], input_record["audit"]
        ).get("ERROR", [])
    )


def _process_analysis_set(
    api: IgvfApi, connection: PConnection, analysis_set: AnalysisSet
) -> tuple[TabularFile | None, set[str], list[str], list[list[str]] | None]:
    """Get annotations file, pseudobulk upload status, upload date, and annotation rows for the provided AnalysisSet record."""
    annotations_file: TabularFile | None = None
    if analysis_set.files is None:
        return None, set(), [], None
    for file_id in analysis_set.files:
        file = utils.retry(num_tries=3)(api.get_by_id)(file_id)
        if file.content_type == "cell annotations":
            annotations_file = cast(TabularFile | None, file.actual_instance)
            if annotations_file is not None and (
                annotations_file.status == "deleted"
                or annotations_file.upload_status == "invalidated"
            ):
                annotations_file = None
            break
    if annotations_file is None or annotations_file.accession is None:
        return None, set(), [], None
    with utils.stream_bytes(api, key=AccessionId(annotations_file.accession)) as f_in:
        metadata_rows = utils.read_tsv_bytes(f_in)
    audit_errors, submission_timestamps = (
        _get_all_pseudobulk_audit_errors_and_timestamps(
            connection=connection,
            analysis_set=analysis_set,
            metadata_rows=metadata_rows,
        )
    )
    return annotations_file, audit_errors, submission_timestamps, metadata_rows


def _get_tracker_row(
    search_result_item: SearchResultItem,
) -> tuple[str, PseudobulkTrackerRow | str, frozenset[str] | None]:
    """For the given AnalysisSet, find required pseudobulk info.

    Args:
        search_result_item: an AnalysisSet that might be pseudobulked.
    Returns:
        analysis_set accession
        PseudobulkTrackerRow with quality info
        required CL_ids for this set.
    """
    # Get analysis set
    analysis_set = cast(AnalysisSet, search_result_item.actual_instance)
    assay = (
        {"NONE"}
        if analysis_set.preferred_assay_titles is None
        else set(analysis_set.preferred_assay_titles)
    )
    if assay & _PREFERRED_ASSAY_TITLES == 0:
        (
            analysis_set.accession,
            f"bad preferred assay titles: {','.join(sorted(assay))}",
            None,
        )

    api = _THREAD_LOCAL_DATA.api
    connection = _THREAD_LOCAL_DATA.connection
    annotations_file, audit_errors, submission_timestamps, metadata_rows = (
        _process_analysis_set(api, connection, analysis_set)
    )

    upload_date = (
        None
        if len(submission_timestamps) == 0
        else max(
            submission_timestamps, key=lambda s: datetime.fromisoformat(s).timestamp()
        )
    )
    num_pseudobulks = len(submission_timestamps)

    # Save annotation file info
    if annotations_file is None:
        return str(analysis_set.accession), "no cell annotations", None
    else:
        try:
            annotations_file_qc = AnnotationsFileQc.qc(
                api=api,
                connection=connection,
                annotations_file_accession=cast(
                    AccessionId | None, annotations_file.accession
                ),
                metadata_rows=metadata_rows,
            )
        except Exception as exception:
            raise ValueError(
                f"Error QC-ing {annotations_file.accession}"
            ) from exception
        pseudobulking_status: str = (
            ",".join(sorted(audit_errors))
            if len(audit_errors) > 0
            else PseudobulkUploadStatus.UNATTEMPTED.value
            if upload_date is None
            else _NO_PSEUDOBULKS
            if num_pseudobulks == 0
            else f"{PseudobulkUploadStatus.COMPLETE.value}: {num_pseudobulks}"
        )
        row: PseudobulkTrackerRow = {
            "principal analysis set accession": f"{analysis_set.accession}",
            "annotation file accession": f"{annotations_file.accession}",
            "lab": f"{annotations_file.lab}",
            "pseudobulking status": pseudobulking_status,
            "processed date": "" if upload_date is None else upload_date,
            "missing annotations columns": ",".join(
                annotations_file_qc.missing_columns
            ),
            "uniform pipeline status": annotations_file_qc.uniform_pipeline_status,
        }
        return str(analysis_set.accession), row, annotations_file_qc.cl_ids


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
    cl_ids_to_update: defaultdict[str, list[str]],
    num_workers: int,
) -> Iterator[PseudobulkTrackerRow]:
    """Iterate over principal analysis sets and yield a PseudobulkTrackerRow with details of their status."""
    if search_result_items is None:
        raise ValueError("No principal analysis sets were found")
    logger.info(f"Scanning analysis sets with {num_workers} workers.")
    lock = Lock()
    with ThreadPoolExecutor(
        max_workers=num_workers, initializer=_init_worker, initargs=(igvf_mode, lock)
    ) as executor:
        for result in executor.map(
            _get_tracker_row, search_result_items, buffersize=2 * num_workers
        ):
            match result:
                case accession, str(reason), _:
                    logger.warning(f"Skipping analysis set {accession} with {reason}.")
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
                        for cl_id in cl_ids:
                            cl_ids_to_update[cl_id].append(accession)
                    yield row


_API: IgvfApi
_CONNECTION: PConnection


def _init_worker(igvf_mode: IgvfMode, lock: Lock):
    """Initialize workers in other processes."""
    _THREAD_LOCAL_DATA.api = utils.open_igvf_api(igvf_mode=igvf_mode)
    _THREAD_LOCAL_DATA.connection = PConnection.new(igvf_mode=igvf_mode, lock=lock)

    # global _API
    # global _CONNECTION
    # _API = utils.open_igvf_api(igvf_mode=igvf_mode)
    # _CONNECTION = PConnection.new(igvf_mode=igvf_mode, lock=lock)
    # utils.fix_igvf_logging(level=logging.WARNING)


def track_pseudobulks(
    *,
    output: Path = Path("pseudobulk_status.tsv"),
    cli_ids_out: Path = Path("missing_cl_ids.txt"),
    compute: bool = True,
    upload: bool = True,
    spreadsheet_url: str = "https://docs.google.com/spreadsheets/d/1iaz5iqB1X_rIE1jnc6y1gKrfQrNDPCrjT-IuVFVadK8/edit",
    google_service_account_json: Path | None = None,
    num_workers: int = _DEFAULT_NUM_WORKERS,
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
    utils.check_access_keys()
    # the root logger is outrageously noisy on info: move it to warning
    logging.getLogger().setLevel(logging.WARNING)
    # get the logger for this tool, and set its level to info
    logger = utils.get_logger_from_file(__file__)
    utils.setup_logger(logger, level=logging.INFO)

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
        cl_ids = defaultdict(list)
        rows_iter = _iter_rows(
            search_results.graph,
            logger=logger,
            igvf_mode=igvf_mode,
            cl_ids_to_update=cl_ids,
            num_workers=_set_num_workers(
                requested_workers=num_workers, default=_DEFAULT_NUM_WORKERS
            ),
        )
        logger.info(f"Writing output to {output}")
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
        missing_ids = sorted(set(cl_ids.keys()) - portal_cl_ids)
        logger.info(f"Writing {len(missing_ids)} missing CL_ids to {cli_ids_out}")
        cli_ids_out.parent.mkdir(exist_ok=True, parents=True)
        with cli_ids_out.open("wt") as f_out:
            f_out.write("missing_id\tanalysis_accessions\n")
            for missing_id in missing_ids:
                analysis_accessions = ",".join(cl_ids[missing_id])
                f_out.write(f"{missing_id}\t{analysis_accessions}\n")

    if upload:
        if (
            google_service_account_json is None
            or not google_service_account_json.exists()
        ):
            raise ValueError(
                "google_service_account_json must point to a valid credentials JSON to upload."
            )
        logger.info(f"Uploading to {spreadsheet_url}")
        upload_dataframe_to_google_sheet(
            output,
            spreadsheet_url=spreadsheet_url,
            worksheet_name="IGVF annotation files",
            credentials_path=google_service_account_json,
        )

    logger.info("Finished.")
