from collections.abc import Iterator
from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

CsrMatrix = scipy.sparse.csc_matrix


def _iter_chunks(size: int, chunk_size: int) -> Iterator[slice[int]]:
    """Yield slices of chunk_size covering range(size); the final slice may be smaller."""
    start = 0
    while start < size:
        end = min(start + chunk_size, size)
        yield slice(start, end)
        start = end


def _describe_obs_and_var(
    adata: ad.AnnData,
    qc_vars: list[str],
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    chunk_size: int = 10000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return pandas matrices of .obs and .var QC columns.

    Both are computed from a single pass over row chunks: slicing the columns of a backed CSR matrix
    reads the whole matrix for every chunk.
    """
    # allocate arrays to store qc results
    num_rows, num_cols = adata.shape
    nonzero = np.empty(shape=(num_rows,), dtype=np.int64)
    # make sum values signed ints just in case there's some weird normalization
    full_sum = np.empty(shape=(num_rows,), dtype=np.int64)
    qc_arrs: list[np.ndarray] = []
    qc_idxs: list[np.ndarray] = []
    for qc_var in qc_vars:
        qc_arrs.append(np.empty(shape=(num_rows,), dtype=np.int64))
        qc_idxs.append(np.flatnonzero(adata.var[qc_var].array))
    # the var sums accumulate across chunks, so accumulate them at full precision
    var_nonzero = np.zeros(shape=(num_cols,), dtype=np.int64)
    var_sum = np.zeros(shape=(num_cols,), dtype=np.float64)

    # fill out obs values chunk-by-chunk, and accumulate var values
    x_sparse = cast(ad.abc.CSRDataset, adata.X)
    for chunk_slice in _iter_chunks(size=num_rows, chunk_size=chunk_size):
        x_chunk = cast(scipy.sparse.csr_matrix, x_sparse[chunk_slice, :])
        nonzero[chunk_slice] = x_chunk.count_nonzero(axis=1)
        full_sum[chunk_slice] = x_chunk.sum(axis=1).flat
        for qc_arr, qc_idx in zip(qc_arrs, qc_idxs, strict=True):
            qc_arr[chunk_slice] = x_chunk[:, qc_idx].sum(axis=1).flat
        var_nonzero += x_chunk.count_nonzero(axis=0)
        var_sum += x_chunk.sum(axis=0, dtype=np.float64).A1

    # finally fill out dataframes and add percentage stats for qc vars
    obs_metrics = pd.DataFrame(index=adata.obs_names)
    obs_metrics[f"n_{var_type}_by_{expr_type}"] = nonzero
    obs_metrics[f"total_{expr_type}"] = full_sum
    with np.errstate(divide="ignore", invalid="ignore"):
        for qc_arr, qc_var in zip(qc_arrs, qc_vars, strict=True):
            obs_metrics[f"total_{expr_type}_{qc_var}"] = qc_arr
            obs_metrics[f"pct_{expr_type}_{qc_var}"] = qc_arr / full_sum * 100

        var_metrics = pd.DataFrame(index=adata.var_names)
        var_metrics[f"n_cells_by_{expr_type}"] = var_nonzero
        var_metrics[f"mean_{expr_type}"] = var_sum / num_rows
        var_metrics[f"pct_dropout_by_{expr_type}"] = (1.0 - var_nonzero / num_rows) * 100.0
        var_metrics[f"total_{expr_type}"] = var_sum.astype(np.int64)

    return obs_metrics, var_metrics


def calculate_qc_metrics(
    adata: ad.AnnData,
    qc_vars: list[str],
    *,
    expr_type: str = "counts",
    var_type: str = "genes",
    chunk_size: int = 10000,
) -> None:
    """Calculate quality control metrics.

    Add QC columns to obs and var.

    This is a replacement of scanpy.pp.calculate_qc_metrics that works with backed adata,
    and with hardcoded inputs: percent_top=None, log1p=False, inplace=True.

    Args:
        adata: AnnData object. Expected to be backed, untested on non-backed adata.
        qc_vars: Keys for boolean columns of .var which identify variables you could want to
            control for (e.g. “ERCC” or “mito”).
        expr_type: Name of kind of values in X.
        var_type: The kind of thing the variables are.
        chunk_size: process (sums, nonzero counts, etc) sparse matrices in chunks this large.
    """
    obs_metrics, var_metrics = _describe_obs_and_var(
        adata=adata,
        qc_vars=qc_vars,
        expr_type=expr_type,
        var_type=var_type,
        chunk_size=chunk_size,
    )
    adata.obs[obs_metrics.columns] = obs_metrics
    adata.var[var_metrics.columns] = var_metrics
