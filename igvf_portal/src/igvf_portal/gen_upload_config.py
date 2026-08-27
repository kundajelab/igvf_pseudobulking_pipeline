import dataclasses
from collections import defaultdict
from collections.abc import (
    Iterator,
    Mapping,
)
from functools import cached_property
from logging import Logger
from pathlib import Path
from types import MappingProxyType
from typing import cast

from igvf_portal import VERSION, utils
from igvf_portal.enums import (
    AnalysisStep,
    ContentType,
    OutputCategory,
)
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.types import (
    AccessionId,
    Alias,
    AnnotationAccessions,
    AnnotationRow,
    CellType,
    IgvfRecord,
    PseudobulkId,
    SampleId,
)


@dataclasses.dataclass(frozen=True, kw_only=True)
class GenUploadConfig:
    """Class for holding config info for gen_upload_script."""

    basedir: Path
    input_file_sets: str | None
    igvf_lookup: IgvfLookup
    compute_md5: bool = True
    dry_run: bool = True
    lab: str = "/labs/anshul-kundaje/"
    award: str = "/awards/HG012069/"
    file_set_type: str = "pseudobulk analysis"
    alias_prefix: str = "anshul-kundaje"
    metadata_path: Path | None = None
    annotations_path: Path | None = None
    logger: Logger
    _content_accessions: dict[tuple[AccessionId, ContentType], set[AccessionId]] = (
        dataclasses.field(default_factory=dict)
    )
    step_1_aliases: defaultdict[Alias, list[Alias]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )
    step_2_aliases: defaultdict[Alias, list[Alias]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )
    step_3_aliases: defaultdict[Alias, list[Alias]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )

    def lookup_record(self, key: Alias | AccessionId) -> IgvfRecord:
        """Convenience method to lookup record within GenUploadConfig."""
        return self.igvf_lookup.lookup_record(key)

    def _lookup_analysis_step_version(self, analysis_step: AnalysisStep) -> list[Alias]:
        """Get aliases for the requested AnalysisStepVersion."""
        step_record = cast(dict[str, object], self.lookup_record(analysis_step.value))
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
                    return self.lookup_record(version_id)["aliases"]
        raise ValueError(
            f"Unable to find version of analysis step {analysis_step} for 'igvf_pseudobulking_pipeline-v{VERSION}'"
        )

    @cached_property
    def input_accessions(self) -> frozenset[AnnotationAccessions] | None:
        """Get AnnotationAccessions objects from input metadata TSV."""
        if self.metadata_path is None:
            return None

        return frozenset(
            AnnotationAccessions.from_csv_row(row)
            for row in utils.iter_csv_rows(
                csv_path=self.metadata_path,
                required_columns=("analysis_set_accession",),
            )
        )

    @cached_property
    def annotations(self) -> Mapping[PseudobulkId, AnnotationRow]:
        """Load pseudobulk annotations TSV and return map from pseudobulk ID to annotations row."""
        # get the required fields
        fields = AnnotationRow.__annotations__
        annotations_path = (
            self.basedir / "cell_name_to_annotation_mapping.tsv"
            if self.annotations_path is None
            else self.annotations_path
        )
        if not annotations_path.exists():
            raise ValueError(f"Annotations path '{annotations_path}' does not exist.")

        def _cast_to_annotations() -> Iterator[AnnotationRow]:
            """Private function to yield the rows of the annotations files as AnnotationRow objects."""
            for row in utils.iter_csv_rows(
                annotations_path, required_columns=AnnotationRow.__required_keys__
            ):
                yield cast(
                    AnnotationRow,
                    {name: fields[name](value) for name, value in row.items()},
                )

        # return lookup dict with keys = pseudobulk ID and values being AnnotationsRows
        annotations = MappingProxyType(
            {
                annotations_row["pseudobulk_id"]: annotations_row
                for annotations_row in _cast_to_annotations()
            }
        )
        self.logger.info(f"Found {len(annotations)} psuedobulk annotations.")
        return annotations

    def _lookup_content_accessions(
        self, analysis_set_accession: AccessionId, content_type: ContentType
    ) -> set[AccessionId]:
        """Query portal for all input files to specified analysis set of specified ContentType.

        If there are no such files, check intermediate analysis sets.
        """
        key = (analysis_set_accession, content_type)
        if key in self._content_accessions:
            return self._content_accessions[key]
        file_accessions = {
            input_file["accession"]
            for input_file in self.lookup_record(analysis_set_accession)["files"]
            if input_file["content_type"] == content_type.value
        }
        if len(file_accessions) == 0:
            # didn't find any references for this analysis set, find references for analysis sets
            # listed in "input_for"
            input_for = self.lookup_record(analysis_set_accession).get(
                "input_for", None
            )
            if input_for is not None:
                file_accessions = {
                    ref
                    for input_for_accession in self.lookup_record(
                        analysis_set_accession
                    )["input_for"]
                    for ref in self._lookup_content_accessions(
                        input_for_accession, content_type=content_type
                    )
                }
        if len(file_accessions) == 0:
            self.logger.warning(
                f"Analysis set {analysis_set_accession} had no files of type {content_type.value}"
            )
        self._content_accessions[key] = file_accessions
        return file_accessions

    def _get_content_type_accession(
        self, input_accession: AnnotationAccessions, content_type: ContentType
    ) -> set[AccessionId]:
        """For a given AnnotationAccessions, get the input files of the requested ContentType."""
        match content_type:
            case ContentType.FRAGMENTS:
                accession_id = input_accession.fragments_file_accession
            case ContentType.MATRIX:
                accession_id = input_accession.matrix_file_accession
            case _:
                accession_id = None
        if accession_id is None:
            # Not explicitly supplied in the input metadata file, look it up from the AnalysisSet
            return self._lookup_content_accessions(
                input_accession.analysis_set_accession,
                content_type=content_type,
            )
        else:
            return {accession_id}

    def _get_all_content_type_accessions(
        self, content_type: ContentType
    ) -> set[AccessionId]:
        """Preferring info in AnnotationAccessions, get all input files of the requested ContentType."""
        if self.input_accessions is None:
            # lookup via file_sets, but this should only happen if metadata is not available, which shouldn't be the case in pipeline
            return {
                fragment_accession
                for input_accession in self.file_sets
                for fragment_accession in self._lookup_content_accessions(
                    input_accession, content_type=content_type
                )
            }
        else:
            # use metadata to get the actual accession info. We should get the right accessions this way.
            return {
                accession
                for input_accession in self.input_accessions
                for accession in self._get_content_type_accession(
                    input_accession=input_accession, content_type=content_type
                )
            }

    @cached_property
    def fragment_accessions(self) -> set[AccessionId]:
        """Get set of unique AccessionIds for input fragments files used by the pipeline."""
        return self._get_all_content_type_accessions(content_type=ContentType.FRAGMENTS)

    @cached_property
    def matrix_accessions(self) -> set[AccessionId]:
        """Get set of unique AccessionIds for input matrix files used by the pipeline."""
        return self._get_all_content_type_accessions(content_type=ContentType.MATRIX)

    def _lookup_aligned_refs(self, accession: AccessionId) -> list[str]:
        """Given the accession ID of an aligned file (e.g. matrix or fragments file) get reference files."""
        return self.lookup_record(accession)["reference_files"]

    @cached_property
    def reference_ids(self) -> set[str]:
        """Get set of reference accession IDs used by aligned input files."""
        return {
            ref_accession
            for input_set in (self.fragment_accessions, self.matrix_accessions)
            for aligned_accession in input_set
            for ref_accession in self._lookup_aligned_refs(aligned_accession)
        }

    @cached_property
    def reference_files(self) -> str:
        """Get comma-separated set of unique reference files used by inputs."""
        return ",".join(self.reference_ids)

    @cached_property
    def assembly(self) -> str:
        """Lookup assembly used in principal analyses."""
        assemblies = {
            self.lookup_record(Alias(reference_file))["assembly"]
            for reference_file in self.reference_ids
        }
        return ",".join(sorted(assemblies))

    @cached_property
    def controlled_access(self) -> bool:
        """Determine if the input set is controlled access."""
        return any(
            self.lookup_record(file_set).get("controlled_access", False)
            for file_set in self.file_sets
        )

    @cached_property
    def annotations_file_accession(self) -> AccessionId | None:
        """Get the accession ID for the annotations file, or return None if it cannot be determined."""
        annotations_files: list[IgvfRecord] = [
            file_record
            for file_set in self.file_sets
            for file_record in self.lookup_record(file_set)["files"]
            if file_record["content_type"] == "cell annotations"
        ]
        if self.metadata_path is not None:
            annotations_md5sum = utils.md5sum(self.metadata_path)
            annotations_files = [
                record
                for record in annotations_files
                if record["md5sum"] == annotations_md5sum
            ]
        record_accessions: set[AccessionId] = {
            record["accession"] for record in annotations_files
        }
        if len(record_accessions) == 1:
            annotations_record = self.lookup_record(record_accessions.pop())
            return annotations_record["accession"]
        else:
            self.logger.warning(
                f"Unable to determine annotations file, got {len(record_accessions)} possible accession IDs."
            )
            return None

    @cached_property
    def annotations_file_alias(self) -> Alias | None:
        """Get a single alias for the annotations file, or return None if it cannot be determined."""
        annotations_accession = self.annotations_file_accession
        if annotations_accession is not None:
            annotations_record = self.lookup_record(annotations_accession)
            return (
                annotations_record["aliases"][0]
                if "aliases" in annotations_record
                else annotations_record["accession"]
            )

    @cached_property
    def analysis_step_versions(self) -> dict[AnalysisStep, list[Alias]]:
        """Get map from AnalysisStep to list of aliases for AnalysisStepVersion"""
        return {
            analysis_step: self._lookup_analysis_step_version(analysis_step)
            for analysis_step in AnalysisStep
        }

    @cached_property
    def input_fragments_aliases(self) -> list[Alias]:
        """Get list of Aliases to fragments files used by the pipeline."""
        return sorted(
            {
                self.lookup_record(fragments_accession)["aliases"][0]
                for fragments_accession in self.fragment_accessions
            }
        )

    @cached_property
    def input_matrices_aliases(self) -> list[Alias]:
        """Get list of Aliases to matrix files used by the pipeline."""
        return sorted(
            {
                self.lookup_record(matrix_accession)["aliases"][0]
                for matrix_accession in self.matrix_accessions
            }
        )

    def derived_from(
        self,
        analysis_step: AnalysisStep,
        output_category: OutputCategory,
        file_set_alias: Alias,
    ) -> list[Alias]:
        """Get list of aliases to files that this data is derived from."""
        annotations_alias_list: list[Alias] = (
            [] if self.annotations_file_alias is None else [self.annotations_file_alias]
        )

        match analysis_step, output_category:
            case AnalysisStep.PSEUDOBULK_ATAC_SEQ, OutputCategory.PSEUDOBULK:
                aliases = sorted(self.input_fragments_aliases + annotations_alias_list)
            case AnalysisStep.PSEUDOBULK_RNA_SEQ, OutputCategory.PSEUDOBULK:
                aliases = sorted(self.input_matrices_aliases + annotations_alias_list)
            case AnalysisStep.PEAK_CALLING, OutputCategory.PSEUDOBULK:
                # peak calling depends on all the ATAC_SEQ pseudobulks in the same pseudobulk file set
                aliases = sorted(set(self.step_1_aliases[file_set_alias]))
            case AnalysisStep.QC, OutputCategory.PSEUDOBULK:
                # pseudobulk QC depends on all the primary pseudobulk files in the same pseudobulk file set
                aliases = sorted(
                    {
                        alias
                        for step_aliases in (
                            self.step_1_aliases,
                            self.step_2_aliases,
                            self.step_3_aliases,
                        )
                        for alias in step_aliases[file_set_alias]
                    }
                )
            case AnalysisStep.QC, _:
                aliases = sorted(
                    {
                        alias
                        for input_aliases in (
                            self.input_fragments_aliases,
                            self.input_matrices_aliases,
                            annotations_alias_list,
                        )
                        if input_aliases is not None
                        for alias in input_aliases
                    }
                )
            case _:
                raise ValueError(
                    f"Invalid combination of AnalysisStep {analysis_step} and OutputCategory {output_category}."
                )
        if aliases is None or len(aliases) == 0:
            raise ValueError(
                f"Unable to find derived_from for analysis_step {analysis_step}"
            )
        return aliases

    @cached_property
    def file_sets(self) -> list[AccessionId]:
        """Return accession IDs of principal analyses for these data."""
        if self.input_file_sets is None or len(self.input_file_sets) == 0:
            return sorted(self._infer_input_file_sets())
        else:
            return cast(list[AccessionId], sorted(self.input_file_sets.split(",")))

    def get_input_file_sets(
        self, cell_name: CellType, subsample: SampleId
    ) -> Iterator[AccessionId]:
        """Return accession IDs of input file_sets for a given pseudobulk"""
        if self.input_accessions is None:
            for file_set in self.file_sets:
                yield file_set
                for input_file_set in self.lookup_record(file_set).get(
                    "input_file_sets", []
                ):
                    if (
                        input_file_set.get("file_set_type", "")
                        == "intermediate analysis"
                    ):
                        yield input_file_set["accession"]
        else:
            for input_accession in self.input_accessions:
                if (
                    input_accession.cell_name == cell_name
                    and input_accession.subsample == subsample
                ):
                    yield input_accession.analysis_set_accession
                    if input_accession.matrix_file_accession is not None:
                        yield self.lookup_record(input_accession.matrix_file_accession)[
                            "file_set"
                        ]["accession"]
                    if input_accession.fragments_file_accession is not None:
                        yield self.lookup_record(
                            input_accession.fragments_file_accession
                        )["file_set"]["accession"]
        if self.annotations_file_accession is not None:
            yield self.lookup_record(self.annotations_file_accession)["file_set"][
                "accession"
            ]

    def _infer_input_file_sets(self) -> set[AccessionId]:
        """Use the IGVF portal to look up parent file sets if they exist, otherwise use provided intermediate"""

        intermediate_path = self.basedir / "analysis_accession_qc_reports"
        intermediate_accessions = (
            {
                AccessionId(_p.name.split("_per_cell_qc.tsv.gz", 1)[0])
                for _p in intermediate_path.glob("*_per_cell_qc.tsv.gz")
            }
            if self.input_accessions is None
            else {
                input_accession.analysis_set_accession
                for input_accession in self.input_accessions
            }
        )
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
        # use "annotation", the cleaned version of cell_name
        return annotations_row["annotation"], annotations_row["subsample"]

    def report_pseudobulk_match(self, pseudobulk_dir: Path) -> None:
        """Report bidirectional match between pseudobulk folders and annotations."""
        folder_ids: set[PseudobulkId] = {
            PseudobulkId(_folder.name)
            for _folder in utils.iter_pseudobulk_dirs(pseudobulk_dir)
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

    def md5sum(self, filepath: Path, chunk_size: int = 2**20) -> str:
        """Compute and md5sum if requested, otherwise return empty string."""
        if self.compute_md5:
            self.logger.info(f"Computing md5 for {filepath}")
            return utils.md5sum(filepath=filepath, chunk_size=chunk_size)
        return ""
