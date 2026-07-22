import abc
import dataclasses
import fnmatch
from pathlib import Path
from typing import cast

import requests

from igvf_portal import utils
from igvf_portal.enums import (
    AnalysisStep,
    OutputCategory,
)
from igvf_portal.gen_upload_config import GenUploadConfig
from igvf_portal.types import (
    Alias,
    PseudobulkId,
    UploadRow,
)


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfUploadBase(abc.ABC):
    """Base class for managing output TSVs for uploading documents to the IGVF portal"""

    output_category: OutputCategory
    match_glob: str | None
    optional: bool = False
    analysis_step: AnalysisStep | None

    @abc.abstractmethod
    def get_row(
        self, check_path: Path, config: GenUploadConfig, doc_aliases: list[Alias]
    ) -> UploadRow | None: ...

    """Function to find the right file using match_glob, update or use doc_aliases, and return the resulting TSV row dict."""

    def matches(self, check_path: Path) -> bool:
        return (
            False
            if self.match_glob is None
            else fnmatch.fnmatch(check_path.name, self.match_glob)
        )

    def get_path(self, folder: Path) -> Path | None:
        """Find the file path that matches match_glob in the supplied folder."""
        if self.match_glob is None:
            raise ValueError("Matching file is not defined for match_glob = None.")
        file_path = next(folder.glob(self.match_glob), None)
        if file_path is None:
            if self.optional:
                return None
            else:
                raise ValueError(f"No output files match '{self.match_glob}'")
        return file_path

    @classmethod
    def _lookup_input_file_sets_alias(cls, config: GenUploadConfig) -> list[Alias]:
        return [
            config.igvf_lookup.lookup_record(file_set)["aliases"][0]
            for file_set in config.file_sets
        ]

    def _get_fileset_alias(self, upload_path: Path, config: GenUploadConfig) -> Alias:
        """Return an alias for the file set that this data is part of."""
        if self.output_category == OutputCategory.PSEUDOBULK:
            # need the pseudobulk file set (for this folder), plus the input file set
            cell_type, sample_id = config.parse_pseudobulk_folder(upload_path)
            # we can only use one file_set, so just pick the first (if there is more than one, it's just aliases anyway)
            file_set = config.file_sets[0]
            return Alias(
                f"{config.alias_prefix}:pseudobulk-{file_set}-{cell_type}-{sample_id}".replace(
                    ",", "_"
                )
            )
        else:
            # we can only use one file_set, so just pick the first (if there is more than one, it's just aliases anyway)
            return self._lookup_input_file_sets_alias(config)[0]

    def _get_file_alias(
        self,
        upload_path: Path,
        config: GenUploadConfig,
        file_set_alias: Alias | None = None,
    ) -> Alias:
        """Return an alias for this data set."""
        if file_set_alias is None:
            file_set_alias = self._get_fileset_alias(
                upload_path=upload_path, config=config
            )
        suffix = upload_path.name.replace(".", "_")
        alias = Alias(f"{file_set_alias}-{suffix}")
        match self.analysis_step:
            case AnalysisStep.PSEUDOBULK_ATAC_SEQ:
                config.step_1_aliases[file_set_alias].append(alias)
            case AnalysisStep.PSEUDOBULK_RNA_SEQ:
                config.step_2_aliases[file_set_alias].append(alias)
            case AnalysisStep.PEAK_CALLING:
                config.step_3_aliases[file_set_alias].append(alias)
        return alias

    def derived_from(
        self, config: GenUploadConfig, file_set_alias: Alias
    ) -> list[Alias]:
        if self.analysis_step is None:
            raise ValueError(
                "Invalid use of derived_from for object with no AnalysisStep"
            )
        return config.derived_from(
            analysis_step=self.analysis_step,
            output_category=self.output_category,
            file_set_alias=file_set_alias,
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfFile(IgvfUploadBase):
    match_glob: str
    file_format: str
    content_type: str
    analysis_step: AnalysisStep
    file_format_specifications: str | None = None

    def _get_row(
        self, check_path: Path, config: GenUploadConfig, **kwargs: str | bool
    ) -> UploadRow | None:
        upload_path = self.get_path(check_path) if check_path.is_dir() else check_path
        if upload_path is None:
            return None
        file_set_alias = self._get_fileset_alias(upload_path=check_path, config=config)
        file_alias = self._get_file_alias(
            upload_path=upload_path, config=config, file_set_alias=file_set_alias
        )

        def _join_aliases(_aliases: list[Alias]) -> str:
            return ",".join(set(_aliases))

        row: UploadRow = {
            "aliases": file_alias,
            "award": config.award,
            "lab": config.lab,
            "derived_manually": False,
            "file_set": file_set_alias,
            "file_format": self.file_format,
            "content_type": self.content_type,
            "md5sum": config.md5sum(upload_path),
            "file_size": upload_path.stat().st_size,
            "submitted_file_name": f"{upload_path.relative_to(config.basedir)}",
            "reference_files": config.reference_files,
            "analysis_step_version": _join_aliases(
                config.analysis_step_versions[self.analysis_step]
            ),
            "derived_from": _join_aliases(
                self.derived_from(config=config, file_set_alias=file_set_alias)
            ),
            **kwargs,
        }
        if self.file_format_specifications is not None:
            row["file_format_specifications"] = self.file_format_specifications

        return row


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class TabularFile(IgvfFile):
    file_format_type: str | None = None

    def get_row(
        self, check_path: Path, config: GenUploadConfig, doc_aliases: list[Alias]
    ) -> UploadRow | None:
        kwargs: dict[str, str | bool] = {"controlled_access": config.controlled_access}
        if self.file_format_type is not None:
            kwargs["file_format_type"] = self.file_format_type
        return super()._get_row(
            check_path=check_path,
            config=config,
            **kwargs,
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class MatrixFile(IgvfFile):
    def get_row(
        self, check_path: Path, config: GenUploadConfig, doc_aliases: list[Alias]
    ) -> UploadRow | None:
        return super()._get_row(
            check_path=check_path,
            config=config,
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class SignalFile(IgvfFile):
    strand_specificity: str

    def get_row(
        self, check_path: Path, config: GenUploadConfig, doc_aliases: list[Alias]
    ) -> UploadRow | None:
        return super()._get_row(
            check_path=check_path,
            config=config,
            strand_specificity=self.strand_specificity,
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfDocument(IgvfUploadBase):
    match_glob: str
    document_type: str
    description: str
    analysis_step: None = None

    def get_row(
        self, check_path: Path, config: GenUploadConfig, doc_aliases: list[Alias]
    ) -> UploadRow | None:
        upload_path = self.get_path(check_path) if check_path.is_dir() else check_path
        if upload_path is None:
            return None
        file_alias = self._get_file_alias(upload_path=upload_path, config=config)
        doc_aliases.append(file_alias)
        cell_type, sample_id = config.parse_pseudobulk_folder(check_path)
        return {
            "aliases": file_alias,
            "award": config.award,
            "lab": config.lab,
            "document_type": self.document_type,
            "description": f"{self.description} for {cell_type} in {sample_id}",
            "attachment": f'{"path": "{upload_path.relative_to(config.basedir)}"}',
        }


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfPseudobulk(IgvfUploadBase):
    output_category: OutputCategory = OutputCategory.PSEUDOBULK
    match_glob: str | None = None
    analysis_step: None = None

    def get_row(
        self,
        check_path: Path,
        config: GenUploadConfig,
        doc_aliases: list[Alias],
    ) -> UploadRow:
        upload_path = check_path if check_path.is_dir() else check_path.parent
        pseudobulk_alias = self._get_fileset_alias(
            upload_path=upload_path, config=config
        )
        # get the annotations for this pseudobulk set in this folder
        pseudobulk_id = PseudobulkId(check_path.name)
        annotations = config.annotations[pseudobulk_id]
        # this format is needed for upload to portal, it corresponds to looking up `term_name` == CL_id
        cl_id = annotations["CL_id"]
        cell_type = f"/sample-terms/{cl_id.replace(':', '_')}/"
        try:
            term_name = config.igvf_lookup.lookup_record(Alias(cell_type))["term_name"]
        except requests.exceptions.HTTPError:
            # note, if this happens, upload will fail. But we can still generate the correct upload script,
            # and manually ask the DACC to add the required SampleTerm
            records = utils.lookup_ontology_by_cl_id(cl_id)
            if records is None:
                raise ValueError(f"Unknown possibly invalid CL_id: {cl_id}")
            term_name = cast(str, records[0]["label"])

        # If cell_qualifier is speicifed in the annotations, use it.
        # Otherwise get cell_qualifier as leftover text after removing term_name from other cell ID
        # description columns.
        cell_qualifier = annotations.get("cell_qualifier", None)
        if cell_qualifier is None:
            cell_qualifier = max(
                (
                    annotations[key].replace(term_name, "").strip()
                    for key in ("cell_name", "cell_description", "CL_term_name")
                ),
                key=lambda _q: len(_q),
            )
            if len(cell_qualifier) == 0:
                # if there is no remainder, just leave cell_qualifier blank
                cell_qualifier = None

        input_file_sets = set(
            config.get_input_file_sets(
                cell_name=annotations["cell_name"],
                subsample=annotations["subsample"],
            )
        )

        row: UploadRow = {
            "aliases": pseudobulk_alias,
            "award": config.award,
            "lab": config.lab,
            "file_set_type": config.file_set_type,
            "cell_type": cell_type,
            "samples": annotations["subsample"],
            "input_file_sets": ",".join(sorted(input_file_sets)),
            "documents": ",".join(sorted(doc_aliases)),
            "merged": False,
        }
        if cell_qualifier is not None:
            row["cell_qualifier"] = cell_qualifier
        return row
