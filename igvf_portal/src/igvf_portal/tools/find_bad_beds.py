import json
import logging
import os
import re
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import local
from typing import cast

from igvf_client.models.pseudobulk_set import PseudobulkSet
from igvf_client.models.search_result_item import SearchResultItem
from igvf_utils.profiles import IgvfSchema

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.constants import VERSION
from igvf_portal.enums import Concurrency, IgvfMode
from igvf_portal.parallel_logger import ParallelLogger
from igvf_portal.register_config import RegisterConfig
from igvf_portal.types import (
    AccessionId,
    Alias,
    IgvfRecord,
    PortalId,
    UploadRow,
)

_BAD_SCORE_PATTERN = re.compile(r"score \([0-9]+\) must be between 0 and 1000")
_BAD_END_PATTERN = re.compile(r"bed->chromEnd\[[0-9]+\] > chromSize\[[0-9]+\]")


_THREAD_LOCAL_DATA: local = local()


def _check_pseudobulk(search_result_item: SearchResultItem) -> None:
    """Check for problems with BED files in this individual pseudobulk."""
    # Get analysis set
    pseudobulk_set = cast(PseudobulkSet, search_result_item.actual_instance)
    pseudobulk_id = PortalId(pseudobulk_set.id)  # ty:ignore[invalid-argument-type]
    connection: PConnection = _THREAD_LOCAL_DATA.connection
    logger = connection.logger
    try:
        pseudobulk_record = connection.lookup_record(pseudobulk_id)
    except Exception:
        logger.warning(f"Unable to look up {pseudobulk_id}, trying by alias")
        aliases = cast(list[Alias], pseudobulk_set.aliases)
        try:
            pseudobulk_record = connection.lookup_record(aliases[0])
        except Exception:
            logger.warning(
                f"Unable to look up {pseudobulk_id} by alias {aliases[0]}, skipping."
            )
            return

    peaks_bed_record: IgvfRecord | None = None
    peaks_big_bed_record: IgvfRecord | None = None
    fragments_tsv_record: IgvfRecord | None = None
    fragments_big_bed_record: IgvfRecord | None = None
    for file_record in pseudobulk_record["files"]:
        if file_record["content_type"] in {
            "peaks",
            "fragments",
        } and utils.file_not_in_portal(file_record):
            logger.warning(
                f"Skipping {file_record['@id']} that is not available in the portal."
            )
            continue
        match file_record["content_type"], file_record["file_format"]:
            case "peaks", "bed":
                peaks_bed_record = connection.lookup_record(file_record["@id"])
            case "peaks", "bigBed":
                peaks_big_bed_record = connection.lookup_record(file_record["@id"])
            case "fragments", "tsv":
                fragments_tsv_record = connection.lookup_record(file_record["@id"])
            case "fragments", "bigBed":
                fragments_big_bed_record = connection.lookup_record(file_record["@id"])
            case _:
                pass

    _check_bed(
        bed_record=peaks_bed_record,
        big_bed_record=peaks_big_bed_record,
        content_type="peaks",
        pseudobulk_id=pseudobulk_id,
        connection=connection,
        logger=logger,
    )

    _check_bed(
        bed_record=fragments_tsv_record,
        big_bed_record=fragments_big_bed_record,
        content_type="fragments",
        pseudobulk_id=pseudobulk_id,
        connection=connection,
        logger=logger,
    )


def _check_bed(
    bed_record: IgvfRecord | None,
    big_bed_record: IgvfRecord | None,
    content_type: str,
    pseudobulk_id: PortalId,
    connection: PConnection,
    logger: ParallelLogger,
) -> None:
    """Check if peaks or fragments file is out-of-spec.

    Checks to insist on:
    -END of bed file to be <= genome length.
    -for peaks files, score <= 1000.
    If corresponding big-bed file does not exist, prep record for creation.
    """
    if bed_record is None:
        logger.info(f"{pseudobulk_id} {content_type} file does not exist.")
        return  # BED file doesn't exist for this pseudobulk, skip check
    try:
        bed_id = bed_record["@id"]
        bed_accession = bed_record["accession"]
        error = cast(
            str,
            json.loads(bed_record.get("validation_error_detail", "{}")).get(
                "validate_files", ""
            ),
        )
    except Exception as exception:
        raise RuntimeError(
            f"Unable to process {pseudobulk_id} {content_type} record:\n"
            f"{json.dumps(bed_record, indent='  ')}"
        ) from exception
    try:
        need_filter = any(
            pattern.search(error) is not None
            for pattern in (_BAD_END_PATTERN, _BAD_SCORE_PATTERN)
        )
        if (not need_filter) and (big_bed_record is not None):
            # don't need to fix this record
            return
        # find the required reference fasta accession
        fasta_accession = _get_fasta_accession(
            bed_record=bed_record, connection=connection, logger=logger
        )
        if fasta_accession is None:
            logger.warning(
                f"Skipping processing {pseudobulk_id} because no valid FASTA could be found."
            )
            return
        # write reference fasta accession
        output: Path = _THREAD_LOCAL_DATA.output
        fasta_accession_file = output / f"{bed_accession}.ref.txt"
        with fasta_accession_file.open("wt") as f_out:
            f_out.write(fasta_accession)

        # write a JSON file with the bed record
        bed_json = output / f"{bed_accession}.bed.json"
        with bed_json.open("wt") as f_out:
            logger.info(f"Dumping bed JSON for {pseudobulk_id} {content_type} file.")
            json.dump(bed_record, f_out)
        # write a JSON file with the preliminary big-bed record
        big_bed_json = output / f"{bed_accession}.bb.json"
        with big_bed_json.open("wt") as f_out:
            logger.info(
                f"Dumping big-bed JSON for {pseudobulk_id} {content_type} file."
            )
            json.dump(
                _bed_record_to_bigbed_record(bed_record)
                if big_bed_record is None
                else big_bed_record,
                f_out,
            )

    except Exception as exception:
        raise RuntimeError(
            f"Processing {pseudobulk_id} {content_type} file '{bed_id}' with record:\n"
            f"{json.dumps(bed_record, indent='  ')}"
        ) from exception


def _get_fasta_accession(
    bed_record: IgvfRecord, connection: PConnection, logger: ParallelLogger
) -> AccessionId | None:
    ref_ids = [
        rec if isinstance(rec, str) else rec["@id"]
        for rec in bed_record["reference_files"]
    ]
    for ref_id in ref_ids:
        record = connection.lookup_record(ref_id, database=True)
        if record["content_type"] == "genome reference":
            if record["status"] == "archived" and record.get("s3_uri", None) is None:
                logger.warning(f"FASTA reference '{ref_id}' is archived")
                return None
            return record["accession"]
    raise ValueError(f"bed_record with id '{bed_record['@id']}' had no FASTA record")


def igvf_record_to_upload_dict(record: IgvfRecord) -> dict[str, object]:
    """Get a dict for uploading to IGVF Portal."""
    return {
        key: val
        for key, val in record.items()
        if key in UploadRow.__annotations__.keys()
    }


def _bed_record_to_bigbed_record(bed_record: IgvfRecord) -> dict[str, object]:
    """Convert the bed record to a corresponding big-bed record."""

    def _to_bb_alias(_alias: str) -> str:
        """Make the big-bed alias distinct but related to the bed alias."""
        if _alias.endswith("_bed_gz") or _alias.endswith("_tsv_gz"):
            return f"{_alias.rsplit('_', 2)[0]}_bb"
        return f"{_alias}_bb"

    bb_dict = {k: v for k, v in bed_record.items()}
    bb_dict["file_format"] = "bigBed"
    bb_dict["aliases"] = sorted(
        {_to_bb_alias(_alias) for _alias in bed_record["aliases"]}
        - set(bed_record["aliases"])
    )
    if bed_record["content_type"] == "peaks":
        bb_dict["file_format_type"] = "bed6+"
    else:
        bb_dict["file_format_type"] = "bed3+"
    match bed_record["submitted_file_name"].rsplit(".", 2):
        case stem, "tsv", "gz":
            record_stem = stem
        case stem, "bed", "gz":
            record_stem = stem
        case stem, something, "gz":
            record_stem = f"{stem}.{something}"
        case _ as name:
            record_stem = name
    bb_dict["submitted_file_name"] = f"{record_stem}.bb"

    return bb_dict


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
    register_config: RegisterConfig,
    output: Path,
    num_workers: int,
) -> None:
    """Iterate over PseudobulkSets and download invalid BED files and their records."""
    if search_result_items is None:
        raise ValueError("No PseudobulkSets were found")

    # Get connection into submit mode:
    connection = register_config.new_connection
    schema: IgvfSchema = connection.profiles.get_profile_from_id("tabular_file")

    output.mkdir(exist_ok=True, parents=True)
    logger.info(f"Scanning PseudobulkSets with {num_workers} workers.")
    with ThreadPoolExecutor(
        max_workers=num_workers,
        initializer=_init_worker,
        initargs=(register_config, schema, output),
    ) as executor:
        for _ in executor.map(
            _check_pseudobulk, search_result_items, buffersize=2 * num_workers
        ):
            pass


def _init_worker(register_config: RegisterConfig, schema: IgvfSchema, output: Path):
    """Initialize workers in other processes."""
    _THREAD_LOCAL_DATA.register_config = register_config
    connection = register_config.new_connection
    connection.profiles.get_profile_from_id(register_config.profile_id)
    connection.profiles._profiles[register_config.cleaned_profile_id] = schema
    _THREAD_LOCAL_DATA.connection = connection
    _THREAD_LOCAL_DATA.output = output


def find_bad_beds(
    *,
    output: Path = Path("."),
    igvf_mode: IgvfMode = IgvfMode.prod,
    continue_on_failed_credentials: bool = False,
    num_workers: int = 16,
    lab: str = "*",
    max_pseudobulks: int | None = None,
) -> None:
    """Find peaks and fragments pseudobulk files that fail audit, and write record JSONs.

    Args:
        output: Path to output folder with data for creating big-beds.
        igvf_mode: Which IGVF server to use: "prod", "staging", or "sandbox"
        continue_on_failed_credentials: If True, when attempting to re-upload a file, if credentials
          cannot be obtained, skip upload and continue. If False, throw exception. Generally this
          results from a file being finalized, and not needing re-upload.
        num_workers: Number of connections to the IGVF Portal. If <= 0, use one per CPU.
        lab: Which lab to limit to when searching for bad pseudobulks. Use '*' for all labs.
        max_pseudobulks: If set to an integer, then quit after checking that number of pseudobulks.
    """
    utils.check_access_keys()
    # the root logger is outrageously noisy on info: move it to warning
    logging.getLogger().setLevel(logging.WARNING)
    # get the logger for this tool, and set its level to info
    logger = utils.get_logger_from_file(__file__)
    utils.setup_logger(logger, level=logging.INFO)
    logger.info(f"Version: {VERSION}")

    register_config = RegisterConfig(
        igvf_mode=igvf_mode,
        profile_id="tabular_file",
        continue_on_failed_credentials=continue_on_failed_credentials,
        concurrency=Concurrency.PROCESS,
    )

    # Log in to IGVF Portal
    api = utils.open_igvf_api(igvf_mode=igvf_mode)

    # Query all pseudobulks with an error
    search_results = api.search(
        type=["PseudobulkSet"],
        limit="all" if max_pseudobulks is None else max_pseudobulks,
        field_filters={
            "status!": "deleted",
            "audit.ERROR.category": "*",
            "lab.@id": lab,
        },
    )
    """
    This would find bad tabular files, but it's not clear how to look for missing big-beds in otherwise okay pseudobulks
    search_results = api.search(
        type=["TabularFile"],
        limit="all",
        field_filters={
            "status!": "deleted",
            "validation_error_detail": "*",
            "lab.@id": "/labs/anshul-kundaje/",
            "file_set.file_set_type": "pseudobulk analysis"
        },
    )
    """

    # Iterate through all queried pseudobulk sets.
    _iter_rows(
        search_results.graph,
        logger=logger,
        register_config=register_config,
        output=output,
        num_workers=_set_num_workers(requested_workers=num_workers, default=12),
    )
