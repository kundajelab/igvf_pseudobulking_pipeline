from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from pseudobulk import sc_qc


@pytest.mark.parametrize(("size", "chunk_size"), [(0, 3), (1, 3), (3, 3), (10, 3), (10, 20)])
def test_iter_chunks_covers_every_index_once(size: int, chunk_size: int) -> None:
    """Chunks must tile range(size) exactly, with no gaps or overlaps."""
    chunks = list(sc_qc._iter_chunks(size=size, chunk_size=chunk_size))
    covered = [idx for chunk in chunks for idx in range(size)[chunk]]
    assert covered == list(range(size))
    assert all(chunk.stop - chunk.start <= chunk_size for chunk in chunks)


def _make_backed_adata(
    tmp_path: Path, num_rows: int, num_cols: int, dtype: type[np.float32 | np.float64] = np.float32
) -> tuple[ad.AnnData, scipy.sparse.csr_matrix]:
    """Make a random backed AnnData, and return it with its count matrix."""
    rng = np.random.default_rng(seed=0)
    x = scipy.sparse.random(
        num_rows,
        num_cols,
        density=0.3,
        format="csr",
        dtype=dtype,
        rng=rng,
        data_rvs=lambda size: rng.integers(1, 20, size=size),
    )
    # include an all-zero cell, to check the percentages of an empty row
    x = scipy.sparse.csr_matrix(scipy.sparse.diags((np.arange(num_rows) > 0).astype(dtype)) @ x)
    x.eliminate_zeros()
    var = pd.DataFrame(
        {
            "mt": np.arange(num_cols) % 5 == 0,
            "ribo": np.arange(num_cols) % 7 == 0,
        },
        index=pd.Index([f"gene{idx}" for idx in range(num_cols)]),
    )
    obs = pd.DataFrame(index=pd.Index([f"cell{idx}" for idx in range(num_rows)]))
    path = tmp_path / "test.h5ad"
    ad.io.write_h5ad(path, ad.AnnData(X=x, obs=obs, var=var))
    return ad.read_h5ad(path, backed="r"), x


# real raw RNA h5ads store counts as float64
@pytest.mark.parametrize("dtype", [np.float64, np.float32])
@pytest.mark.parametrize("chunk_size", [1, 7, 10, 1000])
def test_calculate_qc_metrics_matches_dense(
    tmp_path: Path, chunk_size: int, dtype: type[np.float32 | np.float64]
) -> None:
    """QC metrics must match a direct dense computation, whatever the chunk size.

    sc_qc.calculate_qc_metrics replaces scanpy.pp.calculate_qc_metrics (with percent_top=None,
    log1p=False), so it must also add the same columns, in the same order. Its results were
    checked against scanpy 1.12.4: bit-for-bit identical with float64 or integer counts.
    """
    num_rows, num_cols = 50, 23
    adata, x = _make_backed_adata(tmp_path, num_rows=num_rows, num_cols=num_cols, dtype=dtype)
    dense = x.toarray()
    mt = adata.var["mt"].to_numpy()
    ribo = adata.var["ribo"].to_numpy()

    sc_qc.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], chunk_size=chunk_size)

    obs, var = adata.obs, adata.var
    # the columns that scanpy adds
    assert list(obs.columns) == [
        "n_genes_by_counts",
        "total_counts",
        "total_counts_mt",
        "pct_counts_mt",
        "total_counts_ribo",
        "pct_counts_ribo",
    ]
    assert list(var.columns) == [
        "mt",
        "ribo",
        "n_cells_by_counts",
        "mean_counts",
        "pct_dropout_by_counts",
        "total_counts",
    ]
    total = dense.sum(axis=1)
    np.testing.assert_array_equal(obs["n_genes_by_counts"], (dense != 0).sum(axis=1))
    np.testing.assert_array_equal(obs["total_counts"], total)
    np.testing.assert_array_equal(obs["total_counts_mt"], dense[:, mt].sum(axis=1))
    np.testing.assert_array_equal(obs["total_counts_ribo"], dense[:, ribo].sum(axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        np.testing.assert_allclose(obs["pct_counts_mt"], dense[:, mt].sum(axis=1) / total * 100)
        np.testing.assert_allclose(obs["pct_counts_ribo"], dense[:, ribo].sum(axis=1) / total * 100)

    nonzero = (dense != 0).sum(axis=0)
    np.testing.assert_array_equal(var["n_cells_by_counts"], nonzero)
    np.testing.assert_array_equal(var["total_counts"], dense.sum(axis=0))
    # float32 counts are only accurate to about 1e-7
    rtol = 1e-6 if dtype == np.float32 else 1e-12
    np.testing.assert_allclose(var["mean_counts"], dense.mean(axis=0, dtype=np.float64), rtol=rtol)
    np.testing.assert_allclose(var["pct_dropout_by_counts"], (1.0 - nonzero / num_rows) * 100.0)
