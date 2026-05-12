from enum import (
    Enum,
    StrEnum,
)

from igvf_portal.types import Alias


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
    """Cell with content type"""

    MATRIX = "cell by gene matrix"
    FRAGMENTS = "fragments"


class IgvfMode(StrEnum):
    """Enum for valid IGVF Portal access modes."""

    prod = "prod"
    staging = "staging"
    sandbox = "sandbox"


class OutputCategory(Enum):
    """Enum to describe which kind of file is being uploaded."""

    PRINCIPAL = "principal"
    """File that pertains to entire principal analysis set."""
    INTERMEDIATE = "intermediate"
    """File that pertains to an intermediate analysis set."""
    PSEUDOBULK = "pseudobulk"
    """File that pertains to an individual pseudobulk set."""
