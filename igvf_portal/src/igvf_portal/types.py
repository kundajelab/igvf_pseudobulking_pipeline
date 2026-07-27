import dataclasses
from collections.abc import Mapping
from typing import Literal, NewType, NotRequired, Protocol, TypedDict, cast

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
    CL_term_name: CellType
    subsample: SampleId
    cell_qualifier: NotRequired[CellType]


class UploadRow(TypedDict):
    aliases: str
    award: str
    lab: str
    file_set: NotRequired[str]
    file_set_type: NotRequired[str]
    file_size: NotRequired[int]
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


IgvfRecord = TypedDict(
    "IgvfRecord",
    {
        "@type": list[str],
        "accession": AccessionId,
        "aliases": list[Alias],
        "input_for": NotRequired[list[AccessionId]],
        "input_file_sets": list["IgvfRecord"],
        "file_set": NotRequired["IgvfRecord"],
        "file_set_type": NotRequired["str"],
        "files": list["IgvfRecord"],
        "content_type": str,
        "controlled_access": NotRequired[bool],
        "s3_uri": str,
        "href": str,
        "submitted_file_name": str,
        "status": str,
        "audit": dict[str, object],
        "reference_files": list[str],
        "assembly": str,
        "term_name": NotRequired[str],
        "md5sum": NotRequired[str],
        "summary": NotRequired[str],
    },
)


@dataclasses.dataclass(slots=True, frozen=True, kw_only=True)
class AnnotationAccessions:
    """Class for holding input accession IDs"""

    analysis_set_accession: AccessionId
    cell_name: CellType
    subsample: SampleId
    matrix_file_accession: AccessionId | None = None
    fragments_file_accession: AccessionId | None = None

    @classmethod
    def from_csv_row(cls, row: dict[str, str]) -> AnnotationAccessions:
        analysis_set_accession = row.get("analysis_set_accession", None)
        if analysis_set_accession is None:
            raise ValueError("Must specify analysis_set_accession")
        cell_name = row.get("cell_name", None)
        if cell_name is None:
            raise ValueError("Must specify cell_name")
        subsample = row.get("subsample", None)
        if subsample is None:
            raise ValueError("Must specify subsample")
        return AnnotationAccessions(
            analysis_set_accession=AccessionId(analysis_set_accession),
            cell_name=CellType(cell_name),
            subsample=SampleId(subsample),
            matrix_file_accession=cast(
                AccessionId | None, row.get("matrix_file_accession", None)
            ),
            fragments_file_accession=cast(
                AccessionId | None, row.get("fragments_file_accession", None)
            ),
        )


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


class GeneInfoRow(TypedDict):
    gene_id: str
    gene_name: str
    mt: bool
    ribo: bool


class TssRow(TypedDict):
    gene: str
    transcript: str
    chro: str
    TSS: int
    strand: Literal["+", "-"]


class HasAnnotations(Protocol):
    __annotations__: Mapping[str, object]


class FromTypedDict(Mapping, HasAnnotations): ...
