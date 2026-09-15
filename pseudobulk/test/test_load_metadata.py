import re
from pathlib import Path

import pandas as pd
import pytest

from pseudobulk.utils import load_metadata


def _write_metadata(tmp_path: Path, barcodes_and_accessions: list[tuple[str, str]]) -> Path:
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(
        [
            {
                "barcode_sample": barcode,
                "cell_name": "T cell",
                "cell_description": "a T cell",
                "CL_id": "CL:T",
                "CL_term_name": "T cell",
                "subsample": "s0",
                "analysis_set_accession": accession,
            }
            for barcode, accession in barcodes_and_accessions
        ]
    ).to_csv(metadata, sep="\t", index=False)
    return metadata


def test_load_metadata(tmp_path: Path) -> None:
    metadata = _write_metadata(tmp_path, [("bc1", "ACC1"), ("bc2", "ACC1"), ("bc3", "ACC2")])
    metadata_df = load_metadata(metadata)
    assert list(metadata_df["barcode_sample"]) == ["bc1", "bc2", "bc3"]
    assert (metadata_df["pseudobulk_id"] == "T_cell-s0").all()


@pytest.mark.parametrize(
    "barcodes_and_accessions",
    [
        pytest.param([("bc1", "ACC1"), ("bc2", "ACC1"), ("bc1", "ACC1")], id="same_accession"),
        pytest.param([("bc1", "ACC1"), ("bc2", "ACC1"), ("bc1", "ACC2")], id="other_accession"),
    ],
)
def test_load_metadata_rejects_repeated_barcodes(
    tmp_path: Path, barcodes_and_accessions: list[tuple[str, str]]
) -> None:
    metadata = _write_metadata(tmp_path, barcodes_and_accessions)
    with pytest.raises(ValueError, match=re.escape("1 barcodes on more than one row, e.g. 'bc1'")):
        load_metadata(metadata)
    # the check applies whenever barcodes are read, not just with the default columns
    with pytest.raises(ValueError, match="more than one row"):
        load_metadata(metadata, wanted_cols=["analysis_set_accession", "barcode_sample"])


def _write_annotations(tmp_path: Path, annotations: list[tuple[str, str]]) -> Path:
    """Write metadata with one cell for each (cell_name, subsample)."""
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(
        [
            {
                "barcode_sample": f"bc{idx}",
                "cell_name": cell_name,
                "cell_description": "a cell",
                "CL_id": "CL:1",
                "CL_term_name": "cell",
                "subsample": subsample,
                "analysis_set_accession": "ACC1",
            }
            for idx, (cell_name, subsample) in enumerate(annotations)
        ]
    ).to_csv(metadata, sep="\t", index=False)
    return metadata


@pytest.mark.parametrize(
    ("annotations", "expected_ids"),
    [
        pytest.param(
            # hyphens are kept in cell names, for backwards-compatibility with already submitted
            # pseudobulks, but not in subsamples, so the last hyphen separates the two
            [("CD4-positive T cell", "s0"), ("B cell", "s-1")],
            ["CD4-positive_T_cell-s0", "B_cell-s_1"],
            id="dashes",
        ),
        # all numbers, which pandas would otherwise parse as numbers
        pytest.param([("1", "10"), ("2", "10")], ["1-10", "2-10"], id="numbers"),
    ],
)
def test_load_metadata_pseudobulk_ids(
    tmp_path: Path, annotations: list[tuple[str, str]], expected_ids: list[str]
) -> None:
    """Pseudobulk IDs are made of the sanitized cell name and subsample, joined by a "-"."""
    metadata_df = load_metadata(_write_annotations(tmp_path, annotations))
    assert list(metadata_df["pseudobulk_id"]) == expected_ids
    assert list(metadata_df["cell_name"]) == [cell_name for cell_name, _ in annotations]


@pytest.mark.parametrize("column", ["cell_name", "subsample"])
@pytest.mark.parametrize("value", ["", "  "], ids=["missing", "whitespace"])
def test_load_metadata_rejects_empty_annotations(tmp_path: Path, column: str, value: str) -> None:
    annotations = [("T cell", "s0"), ("B cell", "s0")]
    cell_name, subsample = annotations[1]
    annotations[1] = (value, subsample) if column == "cell_name" else (cell_name, value)
    metadata = _write_annotations(tmp_path, annotations)
    match = f"'{column}' column must not be empty, but is empty on 1 rows"
    with pytest.raises(ValueError, match=re.escape(match)):
        load_metadata(metadata)
    # the check applies when the column is only read to derive the pseudobulk IDs
    with pytest.raises(ValueError, match=re.escape(match)):
        load_metadata(metadata, wanted_cols=["barcode_sample", "pseudobulk_id"])


def test_load_metadata_without_cells(tmp_path: Path) -> None:
    """Metadata with a header but no cells loads, with all of its columns."""
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "barcode_sample\tcell_name\tcell_description\tCL_id\tCL_term_name\tsubsample"
        "\tanalysis_set_accession\n"
    )
    metadata_df = load_metadata(metadata)
    assert len(metadata_df) == 0
    assert {"pseudobulk_id", "cleaned_cell_name", "cleaned_subsample"} <= set(metadata_df.columns)
