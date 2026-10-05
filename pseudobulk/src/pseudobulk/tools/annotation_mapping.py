import logging
from pathlib import Path
from typing import Final

from pseudobulk.utils import load_metadata

_LOGGER_NAME: Final[str] = "annotation-mapping"

ANNOTATION_MAPPING_COLUMNS: Final[tuple[str, ...]] = (
    "pseudobulk_id",
    "cell_name",
    "cleaned_cell_name",
    "CL_id",
    "cell_description",
    "CL_term_name",
    "subsample",
    "cleaned_subsample",
    # NOTE: cell_qualifier is optional, so it is only written when the metadata has it
    "cell_qualifier",
)
"""Columns of the cell name to annotation mapping, in the order they are written."""


def annotation_mapping(*, metadata_loc: Path, output: Path) -> None:
    """Write the mapping from cell name to annotation for every pseudobulk in the metadata.

    The mapping only depends on the metadata, so it is made for every run, including those with no
    RNA data.

    Args:
        metadata_loc: Input annotations metadata file path.
        output: Path to write the mapping TSV to.
    """
    logger = logging.getLogger(name=_LOGGER_NAME)
    cell_name_to_annotation_df = load_metadata(
        metadata_loc, wanted_cols=ANNOTATION_MAPPING_COLUMNS
    ).drop_duplicates()
    cell_name_to_annotation_df.to_csv(output, sep="\t", index=False)
    logger.info(f"Wrote {len(cell_name_to_annotation_df)} annotations to {output}")
