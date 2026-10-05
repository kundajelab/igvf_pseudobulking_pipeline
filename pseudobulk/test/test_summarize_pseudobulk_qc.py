import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from pseudobulk import utils
from pseudobulk.tools.summarize_pseudobulk_qc import summarize_pseudobulk_qc
from pseudobulk.types import FailureAction, FailureHandler

_PSEUDOBULK = "T_cell-s0"
_TSS_ROW_LEN = 4001


def _write_metadata(tmp_path: Path) -> Path:
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(
        [
            {
                "barcode_sample": barcode,
                "cell_name": cell_name,
                "cell_description": f"a {cell_name}",
                "CL_id": f"CL:{cell_name}",
                "CL_term_name": cell_name,
                "subsample": "s0",
                "analysis_set_accession": "ACC1",
            }
            for barcode, cell_name in [("bc1", "T cell"), ("bc2", "T cell"), ("bc3", "B cell")]
        ]
    ).to_csv(metadata, sep="\t", index=False)
    return metadata


def _write_atac_qc(atac_qc_dir: Path) -> None:
    """Write ATAC QC for the two cells of the pseudobulk, a cell of another one, and a non-cell."""
    atac_qc_dir.mkdir()
    rows = [
        # barcode, pseudobulk, num_frags, num_reads, mono-nucleosomal, nucleosome-free
        ("bc1", _PSEUDOBULK, 10, 20, 3, 5),
        ("bc2", _PSEUDOBULK, 5, 5, 1, 2),
        ("bc3", "B_cell-s0", 1000, 2000, 100, 100),
        ("bc4", "null", 1000, 2000, 100, 100),
    ]
    pd.DataFrame(
        [
            {
                "analysis_set_accession": "ACC1",
                "barcode_sample": barcode,
                "pseudobulk_id": pseudobulk,
                "annotated": pseudobulk != "null",
                "num_frags": num_frags,
                "pct_duplicated_reads": (num_reads - num_frags) / num_reads * 100,
                "nucleosomal_signal": (1 + mono) / (1 + nfr),
                "tss_enrichment": 1.0,
                "raw-num_reads": num_reads,
                "raw-num_dup_reads": num_reads - num_frags,
                "raw-mono_nucleosomal_frags": mono,
                "raw-nucleosome_free_frags": nfr,
            }
            for barcode, pseudobulk, num_frags, num_reads, mono, nfr in rows
        ],
        columns=pd.Index(utils.ATAC_QC_COLUMNS),
    ).to_csv(atac_qc_dir / "ACC1.tsv", sep="\t", index=False)
    tss = np.zeros((len(rows), _TSS_ROW_LEN), dtype=np.uint16)
    # bc1 has a center of 11 and flanks of 200, so the pseudobulk has TSS enrichment 1 / (1 + 0.1)
    tss[0, 1995:2006] = 1
    tss[0, 0] = 200
    # the other rows must not be counted
    tss[2:, 2000] = 1000
    scipy.sparse.save_npz(atac_qc_dir / "ACC1_tss_matrix.npz", scipy.sparse.csr_array(tss))


def _write_rna_qc(tmp_path: Path, rna_only_cell: bool = False) -> tuple[Path, Path]:
    """Write RNA QC for the two cells of the pseudobulk, and optionally a cell without ATAC QC."""
    rna_qc = tmp_path / f"{_PSEUDOBULK}.pseudobulked_cell_QC_metrics.tsv"
    barcodes = ["bc1", "bc2", "bc5"] if rna_only_cell else ["bc1", "bc2"]
    pd.DataFrame(
        {
            "analysis_set_accession": "ACC1",
            "barcode_sample": barcodes,
            "annotated": True,
            "pseudobulk_id": _PSEUDOBULK,
            "rna_read_count": [70, 30, 1][: len(barcodes)],
            "gene_count": [3, 2, 1][: len(barcodes)],
            "pct_mito": 10.0,
            "pct_ribo": [20.0, 50.0, 0.0][: len(barcodes)],
        },
        columns=pd.Index(utils.RNA_QC_COLUMNS),
    ).to_csv(rna_qc, sep="\t", index=False)
    counts = tmp_path / f"{_PSEUDOBULK}.pseudobulk_expression.tsv.gz"
    pd.DataFrame(
        {
            "gene_symbol": ["MT-A", "RPL1", "G2", "G3"],
            "mt": [True, False, False, False],
            "ribo": [False, True, False, False],
            "counts": [10, 30, 60, 0],
        },
        index=pd.Index(["E0", "E1", "E2", "E3"], name="gene_id"),
    ).to_csv(counts, sep="\t")
    return rna_qc, counts


def _write_frip(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Write (frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell)."""
    frip_per_cell = tmp_path / f"{_PSEUDOBULK}.frip_per_cell.tsv"
    frip_per_cell.write_text("bc1\t0.4\nbc2\t0\n")
    fragments_per_cell = tmp_path / f"{_PSEUDOBULK}.fragments_per_cell.tsv"
    fragments_per_cell.write_text("bc1\t10\nbc2\t5\n")
    # as written by CALL_PEAKS, cells with no fragments in peaks are missing
    fragments_in_peaks_per_cell = tmp_path / f"{_PSEUDOBULK}.fragments_in_peaks_per_cell.tsv"
    fragments_in_peaks_per_cell.write_text("bc1\t4\n")
    return frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell


def _read_outputs(tmp_path: Path) -> tuple[pd.DataFrame, pd.Series]:
    per_cell = pd.read_csv(tmp_path / "per_cell_qc.tsv.gz", sep="\t")
    summary = pd.read_csv(tmp_path / "pseudobulk_qc.tsv", sep="\t")
    assert len(summary) == 1
    return per_cell, summary.iloc[0]


def test_summarize_pseudobulk_qc(tmp_path: Path) -> None:
    metadata = _write_metadata(tmp_path)
    atac_qc_dir = tmp_path / "atac_qc"
    _write_atac_qc(atac_qc_dir)
    # QC of an analysis set the pseudobulk has no cells in, which is invalid (and has no TSS matrix)
    # so that reading it fails: only the QC of the pseudobulk's own analysis sets should be read
    (atac_qc_dir / "ACC9.tsv").write_text("not\tATAC\nQC\tcolumns\n")
    rna_qc, counts = _write_rna_qc(tmp_path)
    frip_per_cell, fragments_per_cell, fragments_in_peaks_per_cell = _write_frip(tmp_path)

    summarize_pseudobulk_qc(
        pseudobulk=_PSEUDOBULK,
        metadata_loc=metadata,
        atac_qc_dir=atac_qc_dir,
        pseudobulk_qc_out=tmp_path / "per_cell_qc.tsv.gz",
        qc_summary_out=tmp_path / "pseudobulk_qc.tsv",
        rna_qc=rna_qc,
        pseudobulk_counts=counts,
        frip_per_cell=frip_per_cell,
        fragments_per_cell=fragments_per_cell,
        fragments_in_peaks_per_cell=fragments_in_peaks_per_cell,
        failure_actions=[FailureAction.exception],
    )

    per_cell, summary = _read_outputs(tmp_path)
    assert list(per_cell.columns) == [
        "analysis_set_accession",
        "barcode_sample",
        "subsample",
        "rna_read_count",
        "gene_count",
        "pct_mito",
        "pct_ribo",
        "num_frags",
        "pct_duplicated_reads",
        "nucleosomal_signal",
        "tss_enrichment",
        "frip",
    ]
    per_cell = per_cell.set_index("barcode_sample")
    assert sorted(per_cell.index) == ["bc1", "bc2"]
    assert (per_cell["subsample"] == "s0").all()
    assert per_cell["rna_read_count"].to_dict() == {"bc1": 70, "bc2": 30}
    assert per_cell["num_frags"].to_dict() == {"bc1": 10, "bc2": 5}
    assert per_cell["frip"].to_dict() == {"bc1": 0.4, "bc2": 0.0}

    assert summary["pseudobulk"] == _PSEUDOBULK
    assert summary["directory_name"] == _PSEUDOBULK
    assert summary["cell_name"] == "T cell"
    assert summary["subsample"] == "s0"
    assert summary["num_cells"] == 2
    assert summary["rna_read_count"] == 100
    assert summary["gene_count"] == 3
    assert summary["pct_mito"] == pytest.approx(10.0)
    assert summary["pct_ribo"] == pytest.approx(30.0)
    assert summary["num_frags"] == 15
    assert summary["pct_duplicated_reads"] == pytest.approx(10 / 25 * 100)
    assert summary["nucleosomal_signal"] == pytest.approx((1 + 4) / (1 + 7))
    assert summary["tss_enrichment"] == pytest.approx(1 / 1.1)
    assert summary["frip"] == pytest.approx(4 / 15)


def test_summarize_pseudobulk_qc_cell_without_atac(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A cell with RNA QC but no ATAC QC is flagged, and still gets its subsample."""
    metadata = _write_metadata(tmp_path)
    atac_qc_dir = tmp_path / "atac_qc"
    _write_atac_qc(atac_qc_dir)
    rna_qc, counts = _write_rna_qc(tmp_path, rna_only_cell=True)

    with caplog.at_level(logging.WARNING):
        summarize_pseudobulk_qc(
            pseudobulk=_PSEUDOBULK,
            metadata_loc=metadata,
            atac_qc_dir=atac_qc_dir,
            pseudobulk_qc_out=tmp_path / "per_cell_qc.tsv.gz",
            qc_summary_out=tmp_path / "pseudobulk_qc.tsv",
            rna_qc=rna_qc,
            pseudobulk_counts=counts,
            failure_actions=[FailureAction.warning],
        )

    assert "cell sets do not match" in caplog.text
    per_cell, summary = _read_outputs(tmp_path)
    per_cell = per_cell.set_index("barcode_sample")
    assert sorted(per_cell.index) == ["bc1", "bc2", "bc5"]
    assert (per_cell["subsample"] == "s0").all()
    assert per_cell["num_frags"].fillna(-1).to_dict() == {"bc1": 10, "bc2": 5, "bc5": -1}
    assert summary["num_cells"] == 3


def test_summarize_pseudobulk_qc_rna_only(tmp_path: Path) -> None:
    """A pseudobulk with no ATAC data (e.g. RNA-only analysis sets) still knows its subsample."""
    metadata = _write_metadata(tmp_path)
    atac_qc_dir = tmp_path / "atac_qc"
    atac_qc_dir.mkdir()
    rna_qc, counts = _write_rna_qc(tmp_path)

    summarize_pseudobulk_qc(
        pseudobulk=_PSEUDOBULK,
        metadata_loc=metadata,
        atac_qc_dir=atac_qc_dir,
        pseudobulk_qc_out=tmp_path / "per_cell_qc.tsv.gz",
        qc_summary_out=tmp_path / "pseudobulk_qc.tsv",
        rna_qc=rna_qc,
        pseudobulk_counts=counts,
        failure_actions=[FailureAction.exception],
    )

    per_cell, summary = _read_outputs(tmp_path)
    assert per_cell["rna_read_count"].tolist() == [70, 30]
    assert summary["rna_read_count"] == 100
    assert per_cell["subsample"].tolist() == ["s0", "s0"]
    assert per_cell["num_frags"].isna().all()
    # with no fragments there is nothing to count, but the ratios of counts are undefined
    assert summary["num_frags"] == 0
    assert np.isnan(summary["pct_duplicated_reads"])
    assert np.isnan(summary["nucleosomal_signal"])
    assert np.isnan(summary["tss_enrichment"])


@pytest.mark.parametrize("action", list(FailureAction))
def test_failure_handler(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, action: FailureAction
) -> None:
    sentinal = tmp_path / "logs" / "sentinal.log"
    handler = FailureHandler(
        actions=[action], sentinal_file_name=sentinal, logger=logging.getLogger("test")
    )
    if action == FailureAction.exception:
        with pytest.raises(RuntimeError, match="oops"):
            handler.handle_failure("oops")
        return
    with caplog.at_level(logging.WARNING):
        handler.handle_failure("oops")
    if action == FailureAction.warning:
        assert "oops" in caplog.text
    else:
        assert sentinal.read_text() == "oops\n"
