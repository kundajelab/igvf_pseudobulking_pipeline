from typing import NewType, TypeAlias
from typing import TypedDict

PseudobulkId = NewType("PseudobulkId", str)

Alias = NewType("Alias", str)
"""New str type exclusively for aliases."""

CellType = NewType("CellType", str)
"""New str type exclusively for cell types."""

SampleId = NewType("SampleId", str)
"""New str type exclusively for sample IDs."""

AccessionId = NewType("AccessionId", str)
"""New str type exclusively for accession IDs."""

UploadRow: TypeAlias = dict[str, str | bool]
""" Type for rows to upload/register data with the Portal."""


class AnnotationRow(TypedDict):
    """Class with mandatory fields in cell-to-annotation mapping TSV."""

    pseudobulk_id: PseudobulkId
    cell_name: CellType
    annotation: CellType
    CL_id: CellType
    cell_description: CellType
    subsample: SampleId


class IgvfRecord(TypedDict):
    accession: AccessionId
    aliases: list[Alias] | None
    input_for: list[AccessionId]
    input_file_sets: list["IgvfRecord"]
    files: list["IgvfRecord"]
    content_type: str
    controlled_access: bool
