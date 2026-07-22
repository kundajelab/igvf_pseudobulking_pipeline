from pathlib import Path

from igvf_portal import VERSION, utils
from igvf_portal.enums import IgvfMode
from igvf_portal.gen_upload_config import GenUploadConfig
from igvf_portal.igvf_lookup import IgvfLookup
from igvf_portal.upload_state import UploadState


def gen_upload_script(
    basedir: Path,
    *,
    input_file_sets: str | None = None,
    metadata_file: Path | None = None,
    annotations_tsv: Path | None = None,
    compute_md5: bool = True,
    lab: str = "/labs/anshul-kundaje/",
    award: str = "/awards/HG012069/",
    file_set_type: str = "pseudobulk analysis",
    alias_prefix: str = "anshul-kundaje",
    igvf_mode: IgvfMode = IgvfMode.prod,
    dry_run: bool = True,
):
    """Generate TSVs for documents, pseudobulk sets (with document links) and tabular, matrix and signal files.

    A bash upload script `upload.sh` is created inside basedir, and a TSVs used by the upload script are created
    in the subfolder `upload_tsvs`.

    Args:
        basedir: Path to output folder of pseudobulk pipeline
        input_file_sets: comma-separated string of input file sets (analysis sets used to produce the pseudobulks)
        annotations_tsv: CSV/TSV with columns "pseudobulk", "pseudobulk_id", "cell_name", "CL_id", and "cell_description"
        metadata_file: TSV with all annotations needed for for pseudobulks
        outdir: Path to write output files. TSVs for upload submission will be generated in a subfolder "upload-tsvs"
        compute_md5: If true, compute md5 hashes for file (non-document) uploads
        lab: value to use for lab in metadata
        award: value to use for award in metadata
        file_set_type: value to use for file_set_type in metadata
        alias_prefix: value to prefix new aliases with
        reference_files: comma-separated list of reference files (e.g. alignment indices. Should not change much.)
        igvf_mode: use "staging" for testing, "prod" for actual uploads.
        dry_run: if True, do NOT modify the IGVF portal. If False, actually upload pseudobulk results.
    """
    utils.check_access_keys()
    logger = utils.get_logger_from_file(__file__)
    logger.info(f"Version: {VERSION}")

    # store options and helpful info in a big GenUploadConfig object
    config = GenUploadConfig(
        basedir=basedir,
        input_file_sets=input_file_sets,
        metadata_path=metadata_file,
        compute_md5=compute_md5,
        lab=lab,
        award=award,
        file_set_type=file_set_type,
        alias_prefix=alias_prefix,
        dry_run=dry_run,
        igvf_lookup=IgvfLookup.new(igvf_mode=igvf_mode),
        annotations_path=annotations_tsv,
        logger=logger,
    )

    # Create the upload state
    upload_state = UploadState(basedir=basedir, config=config)
    # Get all the row data and commands for uploading to the IGVF portal
    upload_state.get_uploads()
    # Write the upload TSVs and upload command script
    upload_state.write_upload_state()
    # Check that required metadata is present in the IGVF Portal
    upload_state.check_required_metadata()
