import csv
import dataclasses
import shutil
import stat
from pathlib import Path
from collections import defaultdict
from collections.abc import (
    Collection,
    Sequence,
)
from typing import Final

from igvf_portal.config import Config
from igvf_portal.enums import (
    AnalysisStep,
    OutputCategory,
)
from igvf_portal.igvf_uploads import (
    IgvfDocument,
    IgvfPseudobulk,
    IgvfUploadBase,
    MatrixFile,
    SignalFile,
    TabularFile,
)
from igvf_portal.types import (
    Alias,
    UploadRow,
)
from igvf_portal.utils import iter_pseudobulk_dirs


PSEUDOBULK_FILE_DEFINITIONS: Final[tuple[IgvfUploadBase, ...]] = (
    TabularFile(
        analysis_step=AnalysisStep.PSEUDOBULK_ATAC_SEQ,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="fragments.tsv.gz",
        optional=True,
        file_format="tsv",
        content_type="fragments",
        file_format_specifications="buenrostro-bernstein:igvf-single-cell-pipeline-fragment-file-specification",
    ),
    TabularFile(
        analysis_step=AnalysisStep.PEAK_CALLING,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="peaks.narrowPeak.gz",
        optional=True,
        file_format="bed",
        file_format_type="bed6+",
        content_type="peaks",
        file_format_specifications="anshul-kundaje:pseudobulk_peaks_file_format_spec",
    ),
    TabularFile(
        analysis_step=AnalysisStep.PSEUDOBULK_RNA_SEQ,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="pseudobulk_expression.tsv.gz",
        optional=True,
        file_format="tsv",
        content_type="gene quantifications",
        file_format_specifications="anshul-kundaje:pseudobulk_gene_quantifications_file_format_spec",
    ),
    TabularFile(
        analysis_step=AnalysisStep.QC,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="per_cell_qc.tsv.gz",
        file_format="tsv",
        content_type="per-cell quality report",
        file_format_specifications="anshul-kundaje:pseudobulk_per_cell_quality_report_file_format_spec",
    ),
    MatrixFile(
        analysis_step=AnalysisStep.PSEUDOBULK_RNA_SEQ,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="rna_counts_mtx.h5ad",
        optional=True,
        file_format="h5ad",
        content_type="cell by gene matrix",
        file_format_specifications="buenrostro-bernstein:igvf-sc-pipeline-matrix-h5-specification",
    ),
    SignalFile(
        analysis_step=AnalysisStep.PEAK_CALLING,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="peaks_minuslog10pval.bw",
        optional=True,
        file_format="bigWig",
        content_type="ATAC-seq signal p-value",
        strand_specificity="unstranded",
    ),
    SignalFile(
        analysis_step=AnalysisStep.PEAK_CALLING,
        output_category=OutputCategory.PSEUDOBULK,
        match_glob="raw_insertions.bw",
        optional=True,
        file_format="bigWig",
        content_type="ATAC-seq signal",
        strand_specificity="unstranded",
    ),
)

INTERMEDIATE_FILE_DEFINITIONS: Final[tuple[IgvfUploadBase, ...]] = (
    TabularFile(
        analysis_step=AnalysisStep.QC,
        output_category=OutputCategory.INTERMEDIATE,
        match_glob="*_per_cell_qc.tsv.gz",
        file_format="tsv",
        content_type="per-barcode quality report",
        file_format_specifications="anshul-kundaje:pseudobulk_per_barcode_quality_report_file_format_spec",
    ),
)

PRINCIPAL_FILE_DEFINITIONS: Final[tuple[IgvfUploadBase, ...]] = (
    TabularFile(
        analysis_step=AnalysisStep.QC,
        output_category=OutputCategory.PRINCIPAL,
        match_glob="pseudobulk_qc.tsv.gz",
        file_format="tsv",
        content_type="aggregated pseudobulk quality report",
        file_format_specifications="anshul-kundaje:pseudobulk_aggregated_pseudobulk_quality_report_file_format_spec",
    ),
)


@dataclasses.dataclass(kw_only=True, slots=True)
class FileUploadRows:
    tabular_rows: list[UploadRow] = dataclasses.field(default_factory=list)
    matrix_rows: list[UploadRow] = dataclasses.field(default_factory=list)
    signal_rows: list[UploadRow] = dataclasses.field(default_factory=list)


@dataclasses.dataclass(kw_only=True, slots=True)
class UploadState:
    basedir: Path
    config: Config
    document_rows: list[UploadRow] = dataclasses.field(default_factory=list)
    pseudobulk_folders_docs: defaultdict[Path, list[Alias]] = dataclasses.field(
        default_factory=lambda: defaultdict(list)
    )
    analysis_step_file_upload_rows: dict[AnalysisStep, FileUploadRows] = (
        dataclasses.field(
            default_factory=lambda: {step: FileUploadRows() for step in AnalysisStep}
        )
    )
    submission_rows: list[str] = dataclasses.field(default_factory=list)

    def __post_init__(self):
        # Check for missing or extra pseudobulks
        self.config.report_pseudobulk_match(self.pseudobulk_dir)
        # Set initial submission_rows
        self.submission_rows.extend(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                'script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )',
                'pushd "$script_dir" &> /dev/null',
                f'dry_run_arg="{"--dry-run" if self.config.dry_run else ""}"',
                f'igvf_mode="{self.config.igvf_lookup.igvf_mode}"',
            ]
        )

    @property
    def pseudobulk_dir(self) -> Path:
        return self.basedir / "pseudobulks"

    @property
    def intermediate_dir(self) -> Path:
        return self.basedir / "analysis_accession_qc_reports"

    @property
    def upload_tsvs_dir(self) -> Path:
        return self.basedir / "upload_tsvs"

    def _update_rows(
        self,
        file_def: IgvfUploadBase,
        check_path: Path,
        doc_aliases: list[Alias],
        analysis_step: AnalysisStep,
    ) -> None:
        if file_def.analysis_step == analysis_step:
            row = file_def.get_row(
                check_path=check_path, config=self.config, doc_aliases=doc_aliases
            )
            if row is not None:
                file_rows = self.analysis_step_file_upload_rows[analysis_step]
                match file_def:
                    case TabularFile():
                        file_rows.tabular_rows.append(row)
                    case MatrixFile():
                        file_rows.matrix_rows.append(row)
                    case SignalFile():
                        file_rows.signal_rows.append(row)
                    case IgvfDocument():
                        self.document_rows.append(row)

    def _get_pseudobulk_uploads(self, analysis_step: AnalysisStep) -> None:
        """For each pseudobulk folder, get all the uploads in that folder and create a pseudobulk upload."""
        for folder in iter_pseudobulk_dirs(self.pseudobulk_dir):
            doc_aliases = self.pseudobulk_folders_docs[
                folder
            ]  # these aliases will be added to the metadata for this pseudobulk
            for file_def in PSEUDOBULK_FILE_DEFINITIONS:
                self._update_rows(
                    file_def=file_def,
                    check_path=folder,
                    doc_aliases=doc_aliases,
                    analysis_step=analysis_step,
                )

    def _get_intermediate_uploads(self, analysis_step: AnalysisStep) -> None:
        """Get rows and commands for intermediate analysis uploads."""
        doc_aliases = []  # these aliases will not be used
        for check_path in self.intermediate_dir.iterdir():
            for file_def in INTERMEDIATE_FILE_DEFINITIONS:
                if file_def.matches(check_path):
                    self._update_rows(
                        file_def=file_def,
                        check_path=check_path,
                        doc_aliases=doc_aliases,
                        analysis_step=analysis_step,
                    )

    def _get_principal_uploads(self, analysis_step: AnalysisStep) -> None:
        """Get rows and commands for principal analysis uploads."""
        doc_aliases = []  # these aliases will not be used
        for check_path in self.basedir.iterdir():
            for file_def in PRINCIPAL_FILE_DEFINITIONS:
                if file_def.matches(check_path):
                    self._update_rows(
                        file_def=file_def,
                        check_path=check_path,
                        doc_aliases=doc_aliases,
                        analysis_step=analysis_step,
                    )

    def get_uploads(self) -> None:
        """Get rows and commands for all uploads."""
        for analysis_step in AnalysisStep:
            self._get_pseudobulk_uploads(analysis_step)
            self._get_intermediate_uploads(analysis_step)
            self._get_principal_uploads(analysis_step)

    def _write_tsv(
        self,
        upload_type: str,
        rows: Sequence[UploadRow],
        analysis_step: AnalysisStep | None,
        keys: Collection[str] | None = None,
    ) -> None:
        """Write rows as TSV without quoting (important for JSON fields like attachment).
        Return upload line for this dataset
        """
        if len(rows) == 0:
            return

        if keys is None:
            keys = {key for row in rows for key in row.keys()}

        outfile_name = f"{upload_type}.{0 if analysis_step is None else analysis_step.step_num}.tsv"
        outfile = self.upload_tsvs_dir / outfile_name
        outfile.parent.mkdir(exist_ok=True)
        with outfile.open("wt") as f_out:
            writer = csv.DictWriter(
                f_out, fieldnames=keys, delimiter="\t", quoting=csv.QUOTE_NONE
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        self.config.logger.info(f"Wrote {len(rows)} rows to {outfile}")

        step_description = (
            "" if analysis_step is None else f" for step {analysis_step.value}"
        )
        self.submission_rows.append(
            f"1>&2 echo Register {upload_type}{step_description}"
        )
        self.submission_rows.append(
            f'iu_register "$dry_run_arg" -m "$igvf_mode" -p {upload_type} -i "{self.upload_tsvs_dir.name}/{outfile_name}"'
        )

    def write_upload_state(self) -> None:
        """Write the TSVs needed to upload to the IGVF portal"""
        if self.upload_tsvs_dir.exists():
            shutil.rmtree(self.upload_tsvs_dir)
        # Write upload TSVs
        # 1. Documents
        if len(self.document_rows) > 0:
            self._write_tsv("document", self.document_rows, analysis_step=None)
        # 2. Pseudobulks
        if len(self.pseudobulk_folders_docs) > 0:
            igvf_pseudobulk = IgvfPseudobulk()
            pseudobulk_rows = [
                igvf_pseudobulk.get_row(
                    check_path=folder,
                    config=self.config,
                    doc_aliases=doc_aliases,
                )
                for folder, doc_aliases in self.pseudobulk_folders_docs.items()
            ]
            self._write_tsv("pseudobulk_set", pseudobulk_rows, analysis_step=None)
        # then matrix, signal and tabular files in order of step number
        for (
            analysis_step,
            file_upload_rows,
        ) in self.analysis_step_file_upload_rows.items():
            # 3. Matrix files
            if len(file_upload_rows.matrix_rows) > 0:
                self._write_tsv(
                    "matrix_file",
                    file_upload_rows.matrix_rows,
                    analysis_step=analysis_step,
                )
            # 4. Signal files.
            if len(file_upload_rows.signal_rows):
                self._write_tsv(
                    "signal_file",
                    file_upload_rows.signal_rows,
                    analysis_step=analysis_step,
                )
            # 5. Tabular files
            if len(file_upload_rows.tabular_rows):
                self._write_tsv(
                    "tabular_file",
                    file_upload_rows.tabular_rows,
                    analysis_step=analysis_step,
                )

        # add a success message to upload script so that it's clear the entire script succeeded
        self.submission_rows.append("1>&2 echo upload succeeded.")

        # 6. Finally write the upload script
        out = self.basedir / "upload.sh"
        with out.open("wt") as f_out:
            f_out.writelines(f"{row.strip()}\n" for row in self.submission_rows)
        out.chmod(out.stat().st_mode | stat.S_IEXEC)
        self.config.logger.info(f"Wrote submission script to {out}")
