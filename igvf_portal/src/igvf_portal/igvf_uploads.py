import abc
import dataclasses
import fnmatch
from pathlib import Path

from igvf_portal.config import Config
from igvf_portal.enums import (
    AnalysisStep,
    OutputCategory,
)
from igvf_portal.types import (
    Alias,
    PseudobulkId,
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
        self, check_path: Path, config: Config, doc_aliases: list[Alias]
    ) -> dict[str, str] | None: ...

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
    def _lookup_input_file_sets_alias(cls, config: Config) -> list[Alias]:
        return [
            alias
            for file_set in config.file_sets
            for alias in config.igvf_lookup.lookup_aliases(file_set)
        ]

    def _get_fileset_aliases(self, upload_path: Path, config: Config) -> list[Alias]:
        """Return an alias for the file set that this data is part of."""
        match self.output_category:
            case OutputCategory.PSEUDOBULK:
                # get the annotations for this pseudobulk set in this folder
                cell_type, sample_id = config.parse_pseudobulk_folder(upload_path)
                return [
                    Alias(
                        f"{config.alias_prefix}:pseudobulk-{file_set.replace(',', '_')}-{cell_type}-{sample_id}"
                    )
                    for file_set in config.file_sets
                ]
            case OutputCategory.PRINCIPAL | OutputCategory.INTERMEDIATE:
                return self._lookup_input_file_sets_alias(config)

    def _get_aliases(
        self,
        upload_path: Path,
        config: Config,
        file_set_aliases: list[Alias] | None = None,
    ) -> list[Alias]:
        """Return an alias for this data set."""
        if file_set_aliases is None:
            file_set_aliases = self._get_fileset_aliases(
                upload_path=upload_path, config=config
            )
        suffix = upload_path.name.replace(".", "_")
        aliases = [
            Alias(f"{file_set_alias}-{suffix}") for file_set_alias in file_set_aliases
        ]
        match self.analysis_step:
            case AnalysisStep.PSEUDOBULK_ATAC_SEQ:
                config.step_1_aliases.extend(aliases)
            case AnalysisStep.PSEUDOBULK_RNA_SEQ:
                config.step_2_aliases.extend(aliases)
            case AnalysisStep.PEAK_CALLING:
                config.step_3_aliases.extend(aliases)
        return aliases

    def derived_from(self, config: Config) -> list[Alias]:
        if self.analysis_step is None:
            raise ValueError(
                "Invalid use of derived_from for object with no AnalysisStep"
            )
        return config.derived_from(self.analysis_step)


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfFile(IgvfUploadBase):
    match_glob: str
    file_format: str
    content_type: str
    analysis_step: AnalysisStep
    file_format_specifications: str | None = None

    def _get_row(
        self, check_path: Path, config: Config, **kwargs: str
    ) -> dict[str, str] | None:
        upload_path = self.get_path(check_path) if check_path.is_dir() else check_path
        if upload_path is None:
            return None
        file_set_aliases = self._get_fileset_aliases(
            upload_path=check_path, config=config
        )
        file_aliases = self._get_aliases(
            upload_path=upload_path, config=config, file_set_aliases=file_set_aliases
        )

        def _join_aliases(_aliases: list[Alias]) -> str:
            return ",".join(set(_aliases))

        row = {
            "aliases": _join_aliases(file_aliases),
            "award": config.award,
            "lab": config.lab,
            "file_set": _join_aliases(file_set_aliases),
            "file_format": self.file_format,
            "content_type": self.content_type,
            # "assembly": ASSEMBLY,
            "md5sum": config.md5sum(upload_path),
            "submitted_file_name": f"{upload_path.relative_to(config.basedir)}",
            "reference_files": config.reference_files,
            "analysis_step_version": _join_aliases(
                config.analysis_step_versions[self.analysis_step]
            ),
            "derived_from": _join_aliases(self.derived_from(config)),
            **kwargs,
        }
        if self.file_format_specifications is not None:
            row["file_format_specifications"] = self.file_format_specifications

        return row


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class TabularFile(IgvfFile):
    file_format_type: str | None = None

    def get_row(
        self, check_path: Path, config: Config, doc_aliases: list[Alias]
    ) -> dict[str, str] | None:
        kwargs = {"controlled_access": "false"}
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
        self, check_path: Path, config: Config, doc_aliases: list[Alias]
    ) -> dict[str, str] | None:
        return super()._get_row(
            check_path=check_path,
            config=config,
        )


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class SignalFile(IgvfFile):
    strand_specificity: str

    def get_row(
        self, check_path: Path, config: Config, doc_aliases: list[Alias]
    ) -> dict[str, str] | None:
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
        self, check_path: Path, config: Config, doc_aliases: list[Alias]
    ) -> dict[str, str] | None:
        upload_path = self.get_path(check_path) if check_path.is_dir() else check_path
        if upload_path is None:
            return None
        file_aliases = self._get_aliases(upload_path=upload_path, config=config)
        doc_aliases.extend(file_aliases)
        cell_type, sample_id = config.parse_pseudobulk_folder(check_path)
        row = {
            "aliases": ",".join(set(file_aliases)),
            "award": config.award,
            "lab": config.lab,
            "document_type": self.document_type,
            "description": f"{self.description} for {cell_type} in {sample_id}",
            "attachment": f'{"path": "{upload_path.relative_to(config.basedir)}"}',
        }
        return row


@dataclasses.dataclass(frozen=True, kw_only=True, slots=True)
class IgvfPseudobulk(IgvfUploadBase):
    output_category: OutputCategory = OutputCategory.PSEUDOBULK
    match_glob: str | None = None
    analysis_step: None = None

    def get_row(
        self,
        check_path: Path,
        config: Config,
        doc_aliases: list[Alias],
    ) -> dict[str, str]:
        upload_path = check_path if check_path.is_dir() else check_path.parent
        pseudobulk_aliases = self._get_fileset_aliases(
            upload_path=upload_path, config=config
        )
        # get the annotations for this pseudobulk set in this folder
        pseudobulk_id = PseudobulkId(check_path.name)
        annotations = config.annotations[pseudobulk_id]
        return {
            "aliases": ",".join(pseudobulk_aliases),
            "award": config.award,
            "lab": config.lab,
            "file_set_type": config.file_set_type,
            "cell_type": annotations[config.cell_type_key],
            "cell_qualifier": annotations[config.cell_qualifier_key],
            "samples": annotations["subsample"],
            "input_file_sets": ",".join(config.file_sets),
            "documents": ",".join(doc_aliases),
        }
