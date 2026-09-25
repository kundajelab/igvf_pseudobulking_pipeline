import csv
import dataclasses
import logging
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Final

from igvf_client import IgvfApi

from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.gen_upload_config import GenUploadConfig
from igvf_portal.types import AccessionId

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

_MISSING_UNIFORM_PIPELINE_STATUS: Final[str] = "UNKNOWN"
_UNREADABLE_ANNOTATIONS_FILE: Final[str] = "UNREADABLE_ANNOTATIONS"
_NO_ASSEMBLIES_STATUS: Final[str] = "NO_ASSEMBLIES"
_MULTIPLE_ASSEMBLIES_STATUS: Final[str] = "MULTIPLE_ASSEMBLIES"


def _collect_uniform_pipeline_status(
    connection: PConnection,
    accession_ids: Iterable[AccessionId],
) -> set[str] | None:
    statuses = set()
    for accession_id in accession_ids:
        try:
            analysis_set = connection.lookup_record(accession_id)
        except ValueError:
            statuses.add(f"UNPARSEABLE-{accession_id}")
            continue

        status = analysis_set.get("uniform_pipeline_status", None)
        match status:
            case None:
                statuses.add(_MISSING_UNIFORM_PIPELINE_STATUS)
            case _:
                statuses.add(status)
    return None if len(statuses) == 0 else statuses


def _get_assemblies(
    connection: PConnection, annotations_rows: list[list[str]]
) -> tuple[str, ...]:
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False
    ) as temp_file:
        csv_writer = csv.writer(temp_file, delimiter="\t")
        csv_writer.writerows(annotations_rows)
        metadata_path = Path(temp_file.name)
    logger = connection.logger
    original_level = logger.getEffectiveLevel()
    try:
        logger.setLevel(logging.ERROR)
        config = GenUploadConfig(
            basedir=Path("."),
            input_file_sets=None,
            connection=connection,
            metadata_path=metadata_path,
            logger=logger,
        )
        try:
            return config.assemblies
        except KeyError:
            # this data set is broken, return no assemblies
            return ()
    finally:
        logger.setLevel(original_level)
        metadata_path.unlink()


"""
from igvf_portal import utils
from igvf_portal.connection import PConnection
from igvf_portal.annotations_file_qc import AnnotationsFileQc
from igvf_portal.types import AccessionId
api = utils.open_igvf_api("prod")
connection = PConnection("prod")
qc = AnnotationsFileQc.qc(api=api, connection=connection, annotations_file_accession=AccessionId("IGVFFI7967DCIZ"))
qc
"""


@dataclasses.dataclass(slots=True, kw_only=True)
class AnnotationsFileQc:
    missing_columns: frozenset[str]
    on_uniformly_processed_data: set[str] | None
    cl_ids: frozenset[str] | None
    assemblies: tuple[str, ...] = ()
    unreadable_annotations_file: bool = False

    @property
    def uniform_pipeline_status(self) -> str:
        return (
            _UNREADABLE_ANNOTATIONS_FILE
            if self.unreadable_annotations_file
            else _MISSING_UNIFORM_PIPELINE_STATUS
            if self.on_uniformly_processed_data is None
            else ",".join(sorted(self.on_uniformly_processed_data))
        )

    @classmethod
    def get_column_accessions(
        cls, tsv_rows: list[list[str]], column_header: str
    ) -> set[AccessionId]:
        column_index = next(
            (idx for idx, column in enumerate(tsv_rows[0]) if column == column_header),
            None,
        )
        if column_index is None:
            return set()
        return {
            AccessionId(a_id)
            for a_id in (row[column_index].strip() for row in tsv_rows[1:])
            if len(a_id) > 0
        }

    @classmethod
    def qc(
        cls,
        api: IgvfApi,
        connection: PConnection,
        annotations_file_accession: AccessionId | None,
        metadata_rows: list[list[str]] | None = None,
    ) -> AnnotationsFileQc:
        if annotations_file_accession is None:
            raise ValueError("Missing accession for annotations file")
        if metadata_rows is None:
            with utils.stream_bytes(api, key=annotations_file_accession) as f_in:
                tsv_rows = utils.read_tsv_bytes(f_in)
        else:
            tsv_rows = metadata_rows
        if len(tsv_rows) == 0:
            return AnnotationsFileQc(
                missing_columns=_REQUIRED_ANNOTATION_COLUMNS,
                on_uniformly_processed_data=None,
                cl_ids=None,
                unreadable_annotations_file=True,
            )
        missing_columns = set(_REQUIRED_ANNOTATION_COLUMNS.difference(tsv_rows[0]))

        if "analysis_set_accession" in missing_columns:
            on_uniformly_processed_data = None
        else:
            on_uniformly_processed_data = _collect_uniform_pipeline_status(
                connection=connection,
                accession_ids=cls.get_column_accessions(
                    tsv_rows=tsv_rows, column_header="analysis_set_accession"
                ),
            )
            num_assemblies = len(
                _get_assemblies(connection=connection, annotations_rows=tsv_rows)
            )
            match num_assemblies, on_uniformly_processed_data:
                case 1, _:
                    pass
                case 0, None:
                    on_uniformly_processed_data = {_NO_ASSEMBLIES_STATUS}
                case 0, _statuses:
                    on_uniformly_processed_data = _statuses | {_NO_ASSEMBLIES_STATUS}
                case _, None:
                    on_uniformly_processed_data = {_MULTIPLE_ASSEMBLIES_STATUS}
                case _, _statuses:
                    on_uniformly_processed_data = _statuses | {
                        _MULTIPLE_ASSEMBLIES_STATUS
                    }
        if "CL_id" in missing_columns:
            cl_ids = None
        else:
            column_index = next(
                idx for idx, column in enumerate(tsv_rows[0]) if column == "CL_id"
            )
            # replace is needed because some people are putting _ in their CL_id instead of ":"
            cl_ids = frozenset(
                row[column_index].replace("_", ":") for row in tsv_rows[1:]
            )
        return AnnotationsFileQc(
            missing_columns=frozenset(missing_columns),
            on_uniformly_processed_data=on_uniformly_processed_data,
            cl_ids=cl_ids,
        )
