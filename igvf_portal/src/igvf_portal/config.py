from collections.abc import (
    Iterator,
    Mapping,
)
from functools import cached_property
from logging import Logger
from pathlib import Path
from typing import (
    cast,
    Literal,
)
import dataclasses
import hashlib

from igvf_portal import VERSION
from igvf_portal.types import (
    AccessionId,
    Alias,
    AnnotationRow,
    CellType,
    IgvfRecord,
    PseudobulkId,
    SampleId,
)
from igvf_portal.enums import (
    AnalysisStep,
    ContentType,
)
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.utils import iter_pseudobulk_dirs


@dataclasses.dataclass(frozen=True, kw_only=True)
class Config:
    """Class for holding config info."""

    basedir: Path
    input_file_sets: str | None
    igvf_lookup: IgvfLookup
    cell_type_key: Literal["cell_name", "annotation", "CL_id", "cell_description"]
    cell_qualifier_key: Literal["cell_name", "annotation", "CL_id", "cell_description"]
    compute_md5: bool = True
    dry_run: bool = True
    lab: str = "/labs/anshul-kundaje/"
    award: str = "/awards/HG012069/"
    file_set_type: str = "pseudobulk analysis"
    alias_prefix: str = "anshul-kundaje"
    assembly: str = "GRCh38"
    # Reference files (genome + gene annotation)
    reference_files: str = "IGVFFI0653VCGH,IGVFFI9573KOZR"
    annotations: Mapping[PseudobulkId, AnnotationRow]
    logger: Logger
    step_1_aliases: list[Alias] = dataclasses.field(default_factory=list)
    step_2_aliases: list[Alias] = dataclasses.field(default_factory=list)
    step_3_aliases: list[Alias] = dataclasses.field(default_factory=list)

    def _lookup_analysis_step_version(self, analysis_step: AnalysisStep) -> list[Alias]:
        step_record = cast(
            dict[str, object], self.igvf_lookup.lookup_record(analysis_step.value)
        )
        for version_dict in cast(
            list[dict[str, object]], step_record["analysis_step_versions"]
        ):
            for software_versions in cast(
                list[dict[str, str]], version_dict["software_versions"]
            ):
                if (
                    software_versions["name"]
                    == f"igvf_pseudobulking_pipeline-v{VERSION}"
                ):
                    version_id = cast(Alias, version_dict["@id"])
                    return self.igvf_lookup.lookup_aliases(version_id)
        raise ValueError(
            f"Unable to find version of analysis step {analysis_step} for 'igvf_pseudobulking_pipeline-v{VERSION}'"
        )

    @cached_property
    def controlled_access(self) -> bool:
        return any(
            self.igvf_lookup.lookup_record(file_set)["controlled_access"]
            for file_set in self.file_sets
        )

    @cached_property
    def analysis_step_versions(self) -> dict[AnalysisStep, list[Alias]]:
        return {
            analysis_step: self._lookup_analysis_step_version(analysis_step)
            for analysis_step in AnalysisStep
        }

    def derived_from(self, analysis_step: AnalysisStep) -> list[Alias]:
        """Get accession that this data is derived from"""
        match analysis_step:
            case AnalysisStep.PSEUDOBULK_ATAC_SEQ:
                aliases = self.lookup_input_file_set_record_aliases(
                    content_type=ContentType.FRAGMENTS
                )
            case AnalysisStep.PSEUDOBULK_RNA_SEQ:
                aliases = self.lookup_input_file_set_record_aliases(
                    content_type=ContentType.MATRIX
                )
            case AnalysisStep.PEAK_CALLING:
                aliases = self.step_1_aliases
            case AnalysisStep.QC:
                aliases = (
                    self.step_1_aliases + self.step_2_aliases + self.step_3_aliases
                )
        if aliases is None:
            raise ValueError(
                f"Unable to find derived_from for analysis_step {analysis_step}"
            )
        return aliases

    @cached_property
    def file_sets(self) -> list[AccessionId]:
        if self.input_file_sets is None or len(self.input_file_sets) == 0:
            return sorted(self._infer_input_file_sets())
        else:
            return cast(list[AccessionId], sorted(self.input_file_sets.split(",")))

    @classmethod
    def _get_record_file(
        cls, record: IgvfRecord, content_type: ContentType
    ) -> IgvfRecord | None:
        """Get the input file for a IGVF record of the requested ContentType.

        This function assumes there is at most one of that type.
        """
        return next(
            (
                _file
                for _file in record["files"]
                if _file["content_type"] == content_type.value
            ),
            None,
        )

    def _iter_intermediate_records(self, record: IgvfRecord) -> Iterator[IgvfRecord]:
        """Iterate over all accession IDs for"""
        for input_file_set in cast(list[dict[str, object]], record["input_file_sets"]):
            if input_file_set.get("file_set_type", "") == "intermediate analysis":
                yield self.igvf_lookup.lookup_record(
                    cast(AccessionId, input_file_set["accession"])
                )

    def _lookup_input_file_aliases(
        self,
        record: IgvfRecord,
        content_type: ContentType,
        check_intermediate: bool = True,
    ) -> list[Alias] | None:
        """Get the aliases for input files from the input file set record."""
        record_file = self._get_record_file(record, content_type=content_type)
        if record_file is None:
            if check_intermediate:
                return [
                    alias
                    for aliases in (
                        self._lookup_input_file_aliases(
                            intermediate_record,
                            content_type=content_type,
                            check_intermediate=False,
                        )
                        for intermediate_record in self._iter_intermediate_records(
                            record
                        )
                    )
                    if aliases is not None
                    for alias in aliases
                ]
            else:
                return None
        else:
            return record_file["aliases"]

    def lookup_input_file_set_record_aliases(
        self, content_type: ContentType
    ) -> list[Alias] | None:
        """Get all the aliases for files of the specified content type in the input file_set."""
        return [
            alias
            for aliases in (
                self._lookup_input_file_aliases(
                    record=self.igvf_lookup.lookup_record(file_set),
                    content_type=content_type,
                )
                for file_set in self.file_sets
            )
            if aliases is not None
            for alias in aliases
        ]

    def _infer_input_file_sets(self) -> set[AccessionId]:
        """Use the IGVF portal to look up parent file sets if they exist, otherwise use provided intermediate"""

        intermediate_path = self.basedir / "analysis_accession_qc_reports"
        intermediate_accessions = {
            AccessionId(_p.name.split("_per_cell_qc.tsv.gz", 1)[0])
            for _p in intermediate_path.glob("*_per_cell_qc.tsv.gz")
        }
        return self.igvf_lookup.infer_principal_accessions(intermediate_accessions)

    def parse_pseudobulk_folder(
        self, pseudobulk_path: Path
    ) -> tuple[CellType, SampleId]:
        """Extract cell_name and subsample from the pseudobulk folder and cell-name-to-annotations TSV."""
        pseudobulk_id = PseudobulkId(
            (
                pseudobulk_path if pseudobulk_path.is_dir() else pseudobulk_path.parent
            ).name
        )
        annotations_row = self.annotations[pseudobulk_id]
        return annotations_row["cell_name"], annotations_row["subsample"]

    def report_pseudobulk_match(self, pseudobulk_dir: Path) -> None:
        """Report bidirectional match between pseudobulk folders and annotations."""
        folder_ids: set[PseudobulkId] = {
            PseudobulkId(_folder.name)
            for _folder in iter_pseudobulk_dirs(pseudobulk_dir)
        }
        annotation_ids: set[PseudobulkId] = set(self.annotations.keys())

        only_in_folders = folder_ids - annotation_ids  # folders with no lookup entry
        only_in_annotations = (
            annotation_ids - folder_ids
        )  # lookup rows never used by a folder

        self.logger.info("── Cell type match report ──")
        self.logger.info(
            f"  Folders: {len(folder_ids)} distinct pseudobulks | "
            f"Annotation: {len(annotation_ids)} entries"
        )

        if len(only_in_folders) > 0:
            self.logger.error(
                f"  ⚠️  In folders but NOT in annotations ({len(only_in_folders)}) "
                f"-> will cause failures: {sorted(only_in_folders)}"
            )
        if len(only_in_annotations) > 0:
            self.logger.warning(
                f"  ⚠️  In annotations but NOT in any folder ({len(only_in_annotations)}) "
                f"-> unused annotation pseudobulk IDs: {sorted(only_in_annotations)}"
            )
        if len(only_in_folders) == len(only_in_annotations) == 0:
            self.logger.info(
                f"  ✓ Perfect match: all {len(folder_ids)} cell types matched "
                f"in both directions."
            )

    def md5sum(self, filepath: Path, chunk_size: int = 4096) -> str:
        """Compute and md5sum if requested, otherwise return empty string."""
        if self.compute_md5:
            self.logger.info(f"Computing md5 for {filepath}")
            hasher = hashlib.md5()
            with filepath.open("rb") as f:
                while chunk := f.read(chunk_size):
                    hasher.update(chunk)
            return hasher.hexdigest()
        return ""
