import logging
from enum import (
    Enum,
    StrEnum,
)

from igvf_portal.types import Alias, IgvfRecord


class AnalysisStep(Enum):
    """Enum to describe which analysis step of pseudobulking pipeline produced the file.

    To find steps in portal, see:
    https://data.igvf.org/search/?type=AnalysisStep&query=IGVF+Pseudobulking+Step&lab.title=Anshul+Kundaje%2C+Stanford
    """

    PSEUDOBULK_ATAC_SEQ = Alias("anshul-kundaje:igvf-atac-pseudobulking-step")
    """Analysis step for ATAC SEQ pseudobulking.
    See https://data.igvf.org/analysis-steps/eaea908d-3e83-4d49-b111-ef7941188551/
    """
    PSEUDOBULK_RNA_SEQ = Alias("anshul-kundaje:igvf-rna-pseudobulking-step")
    """Analysis step for RNA SEQ pseudobulking.
    See https://data.igvf.org/analysis-steps/cfd4c409-41dc-47f5-bc52-250080e6e3c2/
    """
    PEAK_CALLING = Alias("anshul-kundaje:igvf-pseudobulking-peak-signal-step")
    """Analysis step for peak calling.
    See https://data.igvf.org/analysis-steps/f6c8b4e7-ce6c-4f87-b13b-4301172a172b/
    """
    QC = Alias("anshul-kundaje:igvf-pseudobulking-qc-step")
    """Analysis step for QC.
    See https://data.igvf.org/analysis-steps/673126f6-7367-4e3d-8a7b-832fcefcbebd/
    """

    @property
    def step_num(self) -> int:
        """Ordinal number of this step within the pseudobulking pipeline."""
        match self:
            case AnalysisStep.PSEUDOBULK_ATAC_SEQ:
                return 1
            case AnalysisStep.PSEUDOBULK_RNA_SEQ:
                return 2
            case AnalysisStep.PEAK_CALLING:
                return 3
            case AnalysisStep.QC:
                return 4


class ContentType(Enum):
    """Enum for select IGVF Portal content_types."""

    MATRIX = "cell by gene matrix"
    FRAGMENTS = "fragments"
    PEAKS = "peaks"
    GENOME_REFERENCE = "genome reference"

    @property
    def extension(self) -> str:
        """File extension used for files of this content type."""
        match self:
            case ContentType.MATRIX:
                return "h5ad"
            case ContentType.FRAGMENTS:
                return "bed.gz"
            case ContentType.PEAKS:
                return "tsv.gz"
            case ContentType.GENOME_REFERENCE:
                return "fasta.gz"

    @property
    def file_formats(self) -> frozenset[str]:
        """Allowed file formats used for desired files of this content type."""
        match self:
            case ContentType.MATRIX:
                return frozenset({"h5ad"})
            case ContentType.FRAGMENTS:
                return frozenset({"tsv", "bed"})
            case ContentType.PEAKS:
                return frozenset({"bed"})
            case ContentType.GENOME_REFERENCE:
                return frozenset({"fasta"})

    @staticmethod
    def is_usable(record: IgvfRecord) -> bool:
        """Get whether a file record is neither deleted nor invalidated, whatever its type."""
        return record.get("status") != "deleted" and record.get("upload_status") != "invalidated"

    def is_wanted(self, record: IgvfRecord) -> bool:
        """Get whether a record is of this content type and wanted."""
        return (
            self.is_usable(record)
            and record.get("content_type") == self.value
            and record.get("file_format") in self.file_formats
        )


class IgvfMode(StrEnum):
    """Enum for valid IGVF Portal access modes."""

    prod = "prod"
    staging = "staging"

    @property
    def url(self) -> str:
        """Base API URL of the IGVF Portal for this access mode."""
        match self:
            case IgvfMode.prod:
                return "https://api.data.igvf.org"
            case IgvfMode.staging:
                return "https://api.staging.igvf.org"


class OutputCategory(Enum):
    """Enum to describe which kind of file is being uploaded."""

    PRINCIPAL = "principal"
    """File that pertains to entire principal analysis set."""
    INTERMEDIATE = "intermediate"
    """File that pertains to an intermediate analysis set."""
    PSEUDOBULK = "pseudobulk"
    """File that pertains to an individual pseudobulk set."""


class PseudobulkUploadStatus(Enum):
    """Enum to describe current status of an AnalysisSet's pseudobulks."""

    UNATTEMPTED = "unattempted"
    COMPLETE = "complete"
    NEEDS_FIX = "needs-fix"
    CANNOT_PROCESS = "cannot process"


class MultipleRecordsAction(Enum):
    """Enum to describe action to take if there are multiple records that could be download."""

    KEEP_FILTERED = "keep-filtered"
    KEEP_UNFILTERED = "keep-unfiltered"
    RAISE = "raise"


class LogLevel(Enum):
    """Enum mapping log level names to `logging` level values."""

    info = logging.INFO
    debug = logging.DEBUG
    warning = logging.WARNING
    error = logging.ERROR
    critical = logging.CRITICAL
    fatal = logging.FATAL


class Concurrency(Enum):
    """Enum to describe how work is parallelized: not at all, with threads, or with processes."""

    NONE = None
    THREAD = "Thread"
    PROCESS = "Process"
