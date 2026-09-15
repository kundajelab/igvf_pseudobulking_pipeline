from pathlib import Path

import pandas as pd
import pytest

from pseudobulk.tools.annotation_mapping import annotation_mapping


def _metadata_row(
    barcode: str, cell_name: str, subsample: str, accession: str, **extra: str
) -> dict[str, str]:
    return {
        "barcode_sample": barcode,
        "cell_name": cell_name,
        "cell_description": f"a {cell_name}",
        "CL_id": f"CL:{cell_name}",
        "CL_term_name": cell_name,
        "subsample": subsample,
        "analysis_set_accession": accession,
        **extra,
    }


@pytest.mark.parametrize("with_qualifier", [False, True])
def test_annotation_mapping(tmp_path: Path, with_qualifier: bool) -> None:
    extra = {"cell_qualifier": "activated"} if with_qualifier else {}
    rows = [
        # several cells of the same pseudobulk, across accessions, give one mapping row
        _metadata_row("ACC1_bc0", "T cell", "s0", "ACC1", **extra),
        _metadata_row("ACC1_bc1", "T cell", "s0", "ACC1", **extra),
        _metadata_row("ACC2_bc0", "T cell", "s0", "ACC2", **extra),
        _metadata_row("ACC1_bc2", "B cell", "s1", "ACC1", **extra),
    ]
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(rows).to_csv(metadata, sep="\t", index=False)
    output = tmp_path / "cell_name_to_annotation_mapping.tsv"

    annotation_mapping(metadata_loc=metadata, output=output)

    mapping = pd.read_csv(output, sep="\t", dtype=str)
    expected_columns = [
        "pseudobulk_id",
        "cell_name",
        "cleaned_cell_name",
        "CL_id",
        "cell_description",
        "CL_term_name",
        "subsample",
        "cleaned_subsample",
    ] + (["cell_qualifier"] if with_qualifier else [])
    assert list(mapping.columns) == expected_columns
    assert sorted(mapping["cell_name"]) == ["B cell", "T cell"]
    assert mapping["pseudobulk_id"].is_unique
