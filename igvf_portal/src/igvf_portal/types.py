from typing import NewType, NotRequired, TypedDict

PseudobulkId = NewType("PseudobulkId", str)

Alias = NewType("Alias", str)
"""New str type exclusively for aliases."""

CellType = NewType("CellType", str)
"""New str type exclusively for cell types."""

SampleId = NewType("SampleId", str)
"""New str type exclusively for sample IDs."""

AccessionId = NewType("AccessionId", str)
"""New str type exclusively for accession IDs."""


class AnnotationRow(TypedDict):
    """Class with mandatory fields in cell-to-annotation mapping TSV."""

    pseudobulk_id: PseudobulkId
    cell_name: CellType
    annotation: CellType
    CL_id: CellType
    cell_description: CellType
    subsample: SampleId


class UploadRow(TypedDict):
    aliases: str
    award: str
    lab: str
    file_set: NotRequired[str]
    file_set_type: NotRequired[str]
    file_format: NotRequired[str]
    content_type: NotRequired[str]
    md5sum: NotRequired[str]
    submitted_file_name: NotRequired[str]
    reference_files: NotRequired[str]
    analysis_step_version: NotRequired[str]
    derived_from: NotRequired[str]
    file_format_specifications: NotRequired[str]
    document_type: NotRequired[str]
    description: NotRequired[str]
    attachment: NotRequired[str]
    cell_type: NotRequired[str]
    cell_qualifier: NotRequired[str]
    samples: NotRequired[str]
    input_file_sets: NotRequired[str]
    documents: NotRequired[str]


class IgvfRecord(TypedDict):
    accession: AccessionId
    aliases: list[Alias]
    input_for: NotRequired[list[AccessionId]]
    input_file_sets: list[IgvfRecord]
    files: list[IgvfRecord]
    content_type: str
    controlled_access: bool
    s3_uri: str
    submitted_file_name: str
    status: str
    audit: dict[str, object]


PseudobulkTrackerRow = TypedDict(
    "PseudobulkTrackerRow",
    {
        "principal analysis set accession": str,
        "annotation file accession": str,
        "lab": str,
        "pseudobulking status": str,
        "processed date": str,
        "missing annotations columns": str,
        "uniform pipeline status": str,
    },
)
