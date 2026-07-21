import logging
import os
from collections.abc import Iterator
from pathlib import Path
from typing import (
    Literal,
    cast,
)

from igvf_portal import VERSION
from igvf_portal.config import Config
from igvf_portal.enums import IgvfMode
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.upload_state import UploadState
from igvf_portal.types import (
    PseudobulkId,
    AnnotationRow,
)
from igvf_portal.utils import iter_csv_rows


def _load_annotations(
    tsv_path: Path, logger: logging.Logger
) -> dict[PseudobulkId, AnnotationRow]:
    """Load pseudobulk annotations TSV and return map from pseudobulk ID to annotations row."""
    # get the required fields
    fields = AnnotationRow.__annotations__

    def _cast_to_annotations() -> Iterator[AnnotationRow]:
        """Private function to yield the rows of the annotations files as AnnotationRow objects."""
        for row in iter_csv_rows(tsv_path, required_columns=fields.keys()):
            yield cast(
                AnnotationRow,
                {name: fields[name](value) for name, value in row.items()},
            )

    # return lookup dict with keys = pseudobulk ID and values being AnnotationsRows
    annotations = {
        annotations_row["pseudobulk_id"]: annotations_row
        for annotations_row in _cast_to_annotations()
    }
    logger.info(f"Found {len(annotations)} psuedobulk annotations.")
    return annotations


def gen_upload_script(
    basedir: Path,
    *,
    input_file_sets: str | None = None,
    annotations_tsv: Path | None = None,
    compute_md5: bool = True,
    lab: str = "/labs/anshul-kundaje/",
    award: str = "/awards/HG012069/",
    file_set_type: str = "pseudobulk analysis",
    alias_prefix: str = "anshul-kundaje",
    assembly: str = "GRCh38",
    reference_files: str = "IGVFFI0653VCGH,IGVFFI9573KOZR",
    igvf_mode: IgvfMode = IgvfMode.staging,
    cell_type_key: Literal[
        "cell_name", "annotation", "CL_id", "cell_description"
    ] = "CL_id",
    cell_qualifier_key: Literal[
        "cell_name", "annotation", "CL_id", "cell_description"
    ] = "cell_name",
    dry_run: bool = True,
):
    """Generate TSVs for documents, pseudobulk sets (with document links) and tabular, matrix and signal files.

    A bash upload script `upload.sh` is created inside basedir, and a TSVs used by the upload script are created
    in the subfolder `upload_tsvs`.

    Args:
        basedir: Path to output folder of pseudobulk pipeline
        input_file_sets: comma-separated string of input file sets (analysis sets used to produce the pseudobulks)
        annotations: CSV/TSV with columns "pseudobulk", "pseudobulk_id", "cell_name", "CL_id", and "cell_description"
        outdir: Path to write output files. TSVs for upload submission will be generated in a subfolder "upload-tsvs"
        compute_md5: If true, compute md5 hashes for file (non-document) uploads
        lab: value to use for lab in metadata
        award: value to use for award in metadata
        file_set_type: value to use for file_set_type in metadata
        alias_prefix: value to prefix new aliases with
        assembly: name of assembly for alignment
        reference_files: comma-separated list of reference files (e.g. alignment indices. Should not change much.)
        igvf_mode: use "staging" for testing, "prod" for actual uploads.
        cell_type_key: column to use for "cell_type" metadata in pseudobulks.
        cell_qualifier_key: column to use for "cell_qualifier" metadata in pseudobulks.
        dry_run: if True, do NOT modify the IGVF portal. If False, actually upload pseudobulk results.
    """
    if "IGVF_API_KEY" not in os.environ or "IGVF_SECRET_KEY" not in os.environ:
        raise ValueError("IGVF_API_KEY and IGVF_SECRET_KEY must be set in environment")
    logger = logging.getLogger(name=f"{__package__} generate-igvf-submission-file")
    logger.info(f"Version: {VERSION}")
    annotations = _load_annotations(
        basedir / "cell_name_to_annotation_mapping.tsv"
        if annotations_tsv is None
        else annotations_tsv,
        logger=logger,
    )
    # store options and helpful info in a big Config object
    config = Config(
        basedir=basedir,
        input_file_sets=input_file_sets,
        compute_md5=compute_md5,
        lab=lab,
        award=award,
        file_set_type=file_set_type,
        alias_prefix=alias_prefix,
        assembly=assembly,
        reference_files=reference_files,
        dry_run=dry_run,
        cell_type_key=cell_type_key,
        cell_qualifier_key=cell_qualifier_key,
        igvf_lookup=IgvfLookup.new(igvf_mode=igvf_mode),
        annotations=annotations,
        logger=logger,
    )

    # Create the upload state
    upload_state = UploadState(basedir=basedir, config=config)
    # Get all the row data and commands for uploading to the IGVF portal
    upload_state.get_uploads()
    # Write the upload TSVs and upload command script
    upload_state.write_upload_state()
