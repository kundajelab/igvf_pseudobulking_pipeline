from collections.abc import Callable, Hashable, Sequence

import pytest

from pseudobulk.utils import ATAC_QC_COLUMNS, _select_atac_cols


@pytest.mark.parametrize(
    ("usecols", "expected"),
    [
        (None, list(ATAC_QC_COLUMNS)),
        (
            lambda col: not f"{col}".startswith("raw-"),
            [col for col in ATAC_QC_COLUMNS if not col.startswith("raw-")],
        ),
        (["barcode_sample", "num_frags"], ["barcode_sample", "num_frags"]),
        ([1, 0], [ATAC_QC_COLUMNS[1], ATAC_QC_COLUMNS[0]]),
    ],
)
def test_select_atac_cols(
    usecols: Callable[[Hashable], bool] | Sequence[str] | Sequence[int] | None,
    expected: list[str],
) -> None:
    assert list(_select_atac_cols(usecols)) == expected
