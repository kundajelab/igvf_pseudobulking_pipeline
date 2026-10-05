import logging

import pandas as pd
import pytest

from pseudobulk import utils


def test_merge_rna_and_atac_qc() -> None:
    """Cells in only one of RNA and ATAC QC keep that one's pseudobulk_id and annotated."""
    rna_qc = pd.DataFrame(
        {
            "analysis_set_accession": ["ACC1", "ACC1", "ACC1"],
            "barcode_sample": ["bc1", "bc3", "bc4"],
            "annotated": [True, True, True],
            "pseudobulk_id": ["T_cell-s0", "T_cell-s0", "T_cell-s0"],
            "rna_read_count": [10, 20, 30],
            "gene_count": [5, 6, 7],
            "pct_mito": [1.0, 1.0, 1.0],
            "pct_ribo": [2.0, 2.0, 2.0],
            "found_in_rna": [True, True, True],
        }
    )
    atac_qc = pd.DataFrame(
        {
            "analysis_set_accession": ["ACC1", "ACC1", "ACC1"],
            "barcode_sample": ["bc1", "bc2", "bc4"],
            "pseudobulk_id": ["T_cell-s0", "T_cell-s0", "T_cell-s0"],
            # bc4 has no fragments on an allowed contig
            "annotated": [True, True, False],
            "num_frags": [100, 200, 1],
            "found_in_atac": [True, True, True],
        }
    )
    combined = utils.merge_rna_and_atac_qc(
        identifier="ACC1", rna_qc=rna_qc, atac_qc=atac_qc, logger=logging.getLogger("test")
    ).set_index("barcode_sample")
    assert sorted(combined.index) == ["bc1", "bc2", "bc3", "bc4"]
    # barcode_sample is the index
    assert list(combined.columns) == [
        col for col in utils.FRONT_QC_COLUMNS if col in combined.columns
    ]
    for col in ["pseudobulk_id", "annotated"]:
        assert {f"{col}_rna", f"{col}_atac"}.isdisjoint(combined.columns)
    found_in_rna = {"bc1": True, "bc2": False, "bc3": True, "bc4": True}
    assert combined["found_in_rna"].to_dict() == found_in_rna
    found_in_atac = {"bc1": True, "bc2": True, "bc3": False, "bc4": True}
    assert combined["found_in_atac"].to_dict() == found_in_atac
    assert (combined["pseudobulk_id"] == "T_cell-s0").all()
    assert combined["annotated"].dtype == bool
    assert combined["annotated"].all()


@pytest.mark.parametrize(
    ("has_rna", "has_atac"), [(True, True), (True, False), (False, True), (False, False)]
)
def test_merge_rna_and_atac_qc_column_order(has_rna: bool, has_atac: bool) -> None:
    """The combined QC has the same columns in the same order, whichever QC is present."""
    rna_qc = pd.DataFrame(
        [["ACC1", "bc1", True, "T_cell-s0", 10, 5, 1.0, 2.0, True]],
        columns=pd.Index([*utils.RNA_QC_COLUMNS, "found_in_rna"]),
    )
    # as in combine_accession_qc, which does not load the raw ATAC QC columns
    atac_columns = [col for col in utils.ATAC_QC_COLUMNS if not col.startswith("raw-")]
    atac_qc = pd.DataFrame(
        [["ACC1", "bc1", "T_cell-s0", True, 100, 1.0, 1.0, 1.0, True]],
        columns=pd.Index([*atac_columns, "found_in_atac"]),
    )
    combined = utils.merge_rna_and_atac_qc(
        identifier="ACC1",
        rna_qc=rna_qc if has_rna else rna_qc.iloc[:0],
        atac_qc=atac_qc if has_atac else atac_qc.iloc[:0],
        logger=logging.getLogger("test"),
        raise_on_empty=False,
    )
    assert list(combined.columns) == [
        "analysis_set_accession",
        "barcode_sample",
        "annotated",
        "found_in_rna",
        "found_in_atac",
        "pseudobulk_id",
        "num_frags",
        "pct_duplicated_reads",
        "nucleosomal_signal",
        "tss_enrichment",
        "rna_read_count",
        "gene_count",
        "pct_mito",
        "pct_ribo",
    ]


@pytest.mark.parametrize(("has_rna", "has_atac"), [(True, True), (True, False), (False, True)])
def test_merge_rna_and_atac_qc_row_order(has_rna: bool, has_atac: bool) -> None:
    """The combined QC has its cells in the same order, whichever QC is present."""
    # cells out of order, as when QC files are loaded in filesystem order
    cells = [("ACC2", "bc1"), ("ACC1", "bc2"), ("ACC1", "bc1")]
    rna_qc = pd.DataFrame(
        [[accession, barcode, True, "T_cell-s0", 10, 5, 1.0, 2.0] for accession, barcode in cells],
        columns=pd.Index(utils.RNA_QC_COLUMNS),
    )
    atac_columns = [col for col in utils.ATAC_QC_COLUMNS if not col.startswith("raw-")]
    atac_qc = pd.DataFrame(
        [
            [accession, barcode, "T_cell-s0", True, 100, 1.0, 1.0, 1.0]
            for accession, barcode in cells
        ],
        columns=pd.Index(atac_columns),
    )
    combined = utils.merge_rna_and_atac_qc(
        identifier="test",
        rna_qc=rna_qc if has_rna else rna_qc.iloc[:0],
        atac_qc=atac_qc if has_atac else atac_qc.iloc[:0],
        logger=logging.getLogger("test"),
    )
    assert list(
        zip(combined["analysis_set_accession"], combined["barcode_sample"], strict=True)
    ) == [
        ("ACC1", "bc1"),
        ("ACC1", "bc2"),
        ("ACC2", "bc1"),
    ]
    assert list(combined.index) == [0, 1, 2]
