import dataclasses
import gzip
from collections.abc import Iterable
from io import StringIO
from typing import Final, cast

from igvf_client import IgvfApi
from igvf_client.models.analysis_set import AnalysisSet

from igvf_portal import utils

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


def _collect_uniform_pipeline_status(
    analysis_sets: Iterable[AnalysisSet],
) -> set[str] | None:
    statuses = set()
    for analysis_set in analysis_sets:
        status = analysis_set.uniform_pipeline_status
        match status:
            case None:
                statuses.add(_MISSING_UNIFORM_PIPELINE_STATUS)
            case _:
                statuses.add(status)
    return statuses


@dataclasses.dataclass(slots=True, kw_only=True)
class AnnotationsFileQc:
    missing_columns: frozenset[str]
    on_uniformly_processed_data: set[str] | None
    cl_ids: frozenset[str] | None
    unreadable_annotations_file: bool = False

    @property
    def uniform_pipeline_status(self) -> str:
        return (
            _UNREADABLE_ANNOTATIONS_FILE
            if self.unreadable_annotations_file
            else _MISSING_UNIFORM_PIPELINE_STATUS
            if self.on_uniformly_processed_data is None
            else ",".join(self.on_uniformly_processed_data)
        )

    def qc(api: IgvfApi, annotations_file_accession: str | None) -> AnnotationsFileQc:
        if annotations_file_accession is None:
            raise ValueError("Missing accession for annotations file")
        annotations_compressed_bytes = api.download(annotations_file_accession)
        annotations_decompressed_str = gzip.decompress(
            annotations_compressed_bytes
        ).decode()
        tsv_rows = utils.read_tsv_bytes(StringIO(annotations_decompressed_str))
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
