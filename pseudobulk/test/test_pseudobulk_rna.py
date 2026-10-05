import logging
from pathlib import Path
from typing import cast

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from pseudobulk.tools.pseudobulk_rna import _get_h5ad_paths, _read_rows, pseudobulk_rna

_NUM_GENES = 12
# more cells than sc_qc's default chunk size, so the QC crosses chunk boundaries
_ACCESSION_SIZES = {"ACC1": 10_050, "ACC2": 40}
# (cell_name, subsample) for each annotated cell, chosen by cell index; None leaves it unannotated
_ANNOTATIONS = {
    "ACC1": lambda idx: (
        None if idx % 4 == 3 else ("T cell" if idx % 4 == 0 else "B cell", f"s{idx % 2}")
    ),
    # ACC2 shares "T cell-s0" with ACC1 and has "NK cell-s0" to itself
    "ACC2": lambda idx: None if idx < 5 else ("T cell" if idx < 20 else "NK cell", "s0"),
}


def _metadata_row(barcode: str, cell_name: str, subsample: str, accession: str) -> dict[str, str]:
    return {
        "barcode_sample": barcode,
        "cell_name": cell_name,
        "cell_description": f"a {cell_name}",
        "CL_id": f"CL:{cell_name}",
        "CL_term_name": cell_name,
        "subsample": subsample,
        "analysis_set_accession": accession,
    }


def _counts(adata: ad.AnnData) -> np.ndarray:
    """Get the dense count matrix of an in-memory AnnData."""
    return cast(scipy.sparse.csr_matrix, adata.X).toarray()


def _write_inputs(tmp_path: Path) -> tuple[Path, Path, Path, dict[str, ad.AnnData]]:
    rng = np.random.default_rng(seed=1)
    genes = [f"ENSG{idx:05d}" for idx in range(_NUM_GENES)]
    gene_info = tmp_path / "gene_info.csv"
    pd.DataFrame(
        {
            "gene_name": [f"GENE{idx}" for idx in range(_NUM_GENES)],
            "mt": [idx < 2 for idx in range(_NUM_GENES)],
            "ribo": [2 <= idx < 4 for idx in range(_NUM_GENES)],
        },
        index=pd.Index(genes, name="gene_id"),
    ).to_csv(gene_info)

    h5ad_dir = tmp_path / "h5ads"
    h5ad_dir.mkdir()
    raw: dict[str, ad.AnnData] = {}
    metadata_rows: list[dict[str, str]] = []
    for accession, num_cells in _ACCESSION_SIZES.items():
        barcodes = [f"{accession}_bc{idx}" for idx in range(num_cells)]
        x = scipy.sparse.random(
            num_cells,
            _NUM_GENES,
            density=0.4,
            format="csr",
            dtype=np.float32,
            rng=rng,
            data_rvs=lambda size: rng.integers(1, 50, size=size),
        )
        # real inputs have sparse layers, and could have obsm, which must be carried through too
        adata = ad.AnnData(
            X=x,
            obs=pd.DataFrame(index=pd.Index(barcodes)),
            var=pd.DataFrame(index=pd.Index(genes)),
            layers={"nascent": x * 2},
        )
        adata.obsm["X_embedding"] = rng.random((num_cells, 2))
        ad.io.write_h5ad(h5ad_dir / f"{accession}.h5ad", adata)
        raw[accession] = adata
        for idx, barcode in enumerate(barcodes):
            annotation = _ANNOTATIONS[accession](idx)
            if annotation is None:
                continue
            cell_name, subsample = annotation
            metadata_rows.append(_metadata_row(barcode, cell_name, subsample, accession))
    # an annotated barcode that is in no h5ad, giving an empty pseudobulk
    metadata_rows.append(_metadata_row("ACC2_missing", "Ghost cell", "s0", "ACC2"))
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(metadata_rows).to_csv(metadata, sep="\t", index=False)
    return h5ad_dir, metadata, gene_info, raw


def _expected_pseudobulks(raw: dict[str, ad.AnnData]) -> dict[str, ad.AnnData]:
    """Directly compute the cells of each pseudobulk, in order of accession, then of the h5ad."""
    barcodes: dict[str, list[str]] = {}
    for accession, adata in sorted(raw.items()):
        for idx, barcode in enumerate(adata.obs_names):
            annotation = _ANNOTATIONS[accession](idx)
            if annotation is not None:
                cell_name, subsample = annotation
                barcodes.setdefault(f"{cell_name.replace(' ', '_')}-{subsample}", []).append(
                    barcode
                )
    barcodes["Ghost_cell-s0"] = []
    all_cells = ad.concat(list(raw.values()), axis=0)
    return {
        pseudobulk_id: all_cells[pseudobulk_barcodes, :].copy()
        for pseudobulk_id, pseudobulk_barcodes in barcodes.items()
    }


@pytest.mark.parametrize("num_workers", [1, 3])
def test_pseudobulk_rna(tmp_path: Path, num_workers: int) -> None:
    h5ad_dir, metadata, gene_info, raw = _write_inputs(tmp_path)
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()

    pseudobulk_rna(
        rna_input=h5ad_dir,
        output_dir=output_dir,
        metadata_loc=metadata,
        gene_info=gene_info,
        num_workers=num_workers,
        temp_dir=temp_dir,
    )

    # the temporary per-pseudobulk files are cleaned up
    assert list(temp_dir.iterdir()) == []

    # per-accession QC covers every cell, annotated or not
    for accession, adata in raw.items():
        qc = pd.read_csv(
            output_dir / "rna_qc_reports" / f"{accession}.scRNA_all_cells_QC_metrics.tsv",
            sep="\t",
        )
        dense = _counts(adata)
        total = dense.sum(axis=1)
        assert list(qc["barcode_sample"]) == list(adata.obs_names)
        np.testing.assert_array_equal(qc["rna_read_count"], total)
        np.testing.assert_array_equal(qc["gene_count"], (dense != 0).sum(axis=1))
        # some random cells have no counts, so their percentages are NaN
        with np.errstate(divide="ignore", invalid="ignore"):
            np.testing.assert_allclose(qc["pct_mito"], dense[:, :2].sum(axis=1) / total * 100)
            np.testing.assert_allclose(qc["pct_ribo"], dense[:, 2:4].sum(axis=1) / total * 100)
        annotated = [_ANNOTATIONS[accession](idx) is not None for idx in range(len(adata))]
        assert list(qc["annotated"]) == annotated

    expected = _expected_pseudobulks(raw)
    pseudobulks_dir = output_dir / "pseudobulks"
    assert sorted(path.name for path in pseudobulks_dir.glob("*.rna_counts_mtx.h5ad")) == sorted(
        f"{pseudobulk_id}.rna_counts_mtx.h5ad" for pseudobulk_id in expected
    )
    for pseudobulk_id, expected_adata in expected.items():
        # the large ACC1 finishes last, but its cells must still come first in shared pseudobulks
        actual = ad.read_h5ad(pseudobulks_dir / f"{pseudobulk_id}.rna_counts_mtx.h5ad")
        assert list(actual.obs_names) == list(expected_adata.obs_names)
        assert list(actual.var_names) == list(expected_adata.var_names)
        np.testing.assert_array_equal(_counts(actual), _counts(expected_adata))
        assert list(actual.layers) == ["nascent"]
        np.testing.assert_array_equal(
            cast(scipy.sparse.csr_matrix, actual.layers["nascent"]).toarray(),
            cast(scipy.sparse.csr_matrix, expected_adata.layers["nascent"]).toarray(),
        )
        assert list(actual.obsm) == ["X_embedding"]
        np.testing.assert_array_equal(
            actual.obsm["X_embedding"], expected_adata.obsm["X_embedding"]
        )
        assert (actual.obs["pseudobulk_id"] == pseudobulk_id).all()

        expression = pd.read_csv(
            pseudobulks_dir / f"{pseudobulk_id}.pseudobulk_expression.tsv.gz",
            sep="\t",
            index_col=0,
        )
        np.testing.assert_array_equal(
            expression["counts"],
            _counts(expected_adata).sum(axis=0),
        )

        pseudobulk_qc = pd.read_csv(
            output_dir / "rna_qc_reports" / f"{pseudobulk_id}.pseudobulked_cell_QC_metrics.tsv",
            sep="\t",
        )
        assert list(pseudobulk_qc["barcode_sample"]) == list(expected_adata.obs_names)


def test_pseudobulk_rna_is_reproducible(tmp_path: Path) -> None:
    """Outputs must be byte-identical whatever the number of workers, and so the order of tasks."""
    h5ad_dir, metadata, gene_info, _ = _write_inputs(tmp_path)
    output_dirs: list[Path] = []
    for num_workers in (1, 3):
        output_dir = tmp_path / f"output_{num_workers}"
        output_dir.mkdir()
        temp_dir = tmp_path / f"temp_{num_workers}"
        temp_dir.mkdir()
        pseudobulk_rna(
            rna_input=h5ad_dir,
            output_dir=output_dir,
            metadata_loc=metadata,
            gene_info=gene_info,
            num_workers=num_workers,
            temp_dir=temp_dir,
        )
        output_dirs.append(output_dir)

    def output_files(output_dir: Path) -> dict[Path, bytes]:
        return {
            path.relative_to(output_dir): path.read_bytes()
            for path in sorted(output_dir.rglob("*"))
            if path.is_file()
        }

    single_worker, multi_worker = (output_files(output_dir) for output_dir in output_dirs)
    assert list(single_worker) == list(multi_worker)
    assert len(single_worker) > 0
    for path, contents in single_worker.items():
        assert contents == multi_worker[path], f"{path} differs between worker counts"


def test_get_h5ad_paths_from_file_of_files(tmp_path: Path) -> None:
    h5ads = [tmp_path / "ACC1.h5ad", tmp_path / "ACC2.h5ad"]
    file_of_files = tmp_path / "h5ads.txt"
    file_of_files.write_text("".join(f"{h5ad}\n" for h5ad in h5ads))
    assert list(_get_h5ad_paths(file_of_files, logger=logging.getLogger("test"))) == h5ads


@pytest.mark.parametrize("row_idx", [[0, 2, 3], []], ids=["rows", "no_rows"])
def test_read_rows(tmp_path: Path, row_idx: list[int]) -> None:
    """Each kind of element in an h5ad has just the wanted rows read back."""
    rng = np.random.default_rng(seed=0)
    x = scipy.sparse.random(5, 4, density=0.5, format="csr", dtype=np.float64, rng=rng)
    adata = ad.AnnData(
        X=x,
        obs=pd.DataFrame(index=pd.Index([f"bc{idx}" for idx in range(5)])),
        layers={"csc": scipy.sparse.csc_matrix(x), "dense": x.toarray()},
    )
    adata.obsm["array"] = rng.random((5, 2))
    frame = pd.DataFrame({"a": range(5), "b": list("vwxyz")}, index=adata.obs_names)
    adata.obsm["frame"] = frame
    path = tmp_path / "test.h5ad"
    ad.io.write_h5ad(path, adata)

    idx = np.array(row_idx, dtype=np.intp)
    with h5py.File(path, "r") as h5ad:
        x_rows = _read_rows(h5ad["X"], idx)
        csc_rows = _read_rows(h5ad["layers/csc"], idx)
        dense_rows = _read_rows(h5ad["layers/dense"], idx)
        array_rows = _read_rows(h5ad["obsm/array"], idx)
        frame_rows = _read_rows(h5ad["obsm/frame"], idx)

    assert isinstance(x_rows, scipy.sparse.csr_matrix | scipy.sparse.csr_array)
    np.testing.assert_array_equal(x_rows.toarray(), x.toarray()[idx])
    assert isinstance(csc_rows, scipy.sparse.csc_matrix | scipy.sparse.csc_array)
    np.testing.assert_array_equal(csc_rows.toarray(), x.toarray()[idx])
    assert isinstance(dense_rows, np.ndarray)
    np.testing.assert_array_equal(dense_rows, x.toarray()[idx])
    assert isinstance(array_rows, np.ndarray)
    np.testing.assert_array_equal(array_rows, adata.obsm["array"][idx])
    assert isinstance(frame_rows, pd.DataFrame)
    pd.testing.assert_frame_equal(frame_rows, frame.iloc[idx])
