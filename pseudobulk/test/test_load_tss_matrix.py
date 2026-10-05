from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from pseudobulk.utils import _load_tss_matrix


def test_load_tss_matrix_filters_rows(tmp_path: Path) -> None:
    tss = np.arange(12, dtype=np.uint16).reshape(3, 4)
    matrix_file = tmp_path / "ACC1_tss_matrix.npz"
    scipy.sparse.save_npz(matrix_file, scipy.sparse.csr_array(tss))

    loaded = _load_tss_matrix(matrix_file, row_is_wanted=None)
    np.testing.assert_array_equal(loaded.toarray(), tss)
    filtered = _load_tss_matrix(matrix_file, row_is_wanted=pd.Series([True, False, True]))
    assert isinstance(filtered, scipy.sparse.csr_array)
    np.testing.assert_array_equal(filtered.toarray(), tss[[0, 2]])


def test_load_tss_matrix_rejects_other_formats(tmp_path: Path) -> None:
    matrix_file = tmp_path / "ACC1_tss_matrix.npz"
    scipy.sparse.save_npz(matrix_file, scipy.sparse.csr_matrix(np.eye(2, dtype=np.uint16)))
    with pytest.raises(TypeError, match="csr_matrix"):
        _load_tss_matrix(matrix_file, row_is_wanted=None)
