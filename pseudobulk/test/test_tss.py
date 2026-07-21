from collections.abc import Iterable

import numpy as np
import pytest

from pseudobulk.types import COUNTS_DTYPE
from pseudobulk.types import COUNTS_ARRAY
from pseudobulk.types import POS_DTYPE
from pseudobulk.tss import bisect_left as bisect_left_cython
from pseudobulk.tss import bisect_right as bisect_right_cython
from pseudobulk.tss import update_tss_insertions as update_tss_insertions_cython


def bisect_left(arr: Iterable[int], val: int) -> int:
    """Wrap bisect_left_cython to convert from easily-definable python types to numpy."""
    return int(bisect_left_cython(np.array(arr, dtype=np.int32), np.int32(val)))


def bisect_right(arr: Iterable[int], val: int) -> int:
    """Wrap bisect_right_cython to convert from easily-definable python types to numpy."""
    return int(bisect_right_cython(np.array(arr, dtype=np.int32), np.int32(val)))


def update_tss_insertions(
    tss_insertions: Iterable[int],
    position: int,
    tss_vec: Iterable[int],
    strand_vec: Iterable[int],
) -> COUNTS_ARRAY:
    np_tss_insertions = np.array(tss_insertions, dtype=COUNTS_DTYPE)
    update_tss_insertions_cython(
        tss_insertions=np_tss_insertions,
        position=POS_DTYPE.type(position),
        tss_vec=np.array(tss_vec, dtype=POS_DTYPE),
        strand_vec=np.array(strand_vec, dtype=POS_DTYPE),
    )
    return np_tss_insertions


def test_bisect_left():
    """Test that home-rolled cython bisect_left is correct."""
    assert bisect_left([], 5) == 0
    assert bisect_left([0, 5, 10], 0) == 0
    assert bisect_left([0, 0, 0, 5, 10], 0) == 0
    assert bisect_left([0, 5, 10], -1) == 0
    assert bisect_left([0, 5, 10], 3) == 1
    assert bisect_left([0, 5, 10], 5) == 1
    assert bisect_left([0, 5, 5, 5, 10], 5) == 1
    assert bisect_left([0, 5, 10], 10) == 2
    assert bisect_left([0, 5, 10, 10, 10], 10) == 2
    assert bisect_left([0, 5, 10], 11) == 3


def test_bisect_right():
    """Test that home-rolled cython bisect_right is correct."""
    assert bisect_right([], 5) == 0
    assert bisect_right([0, 5, 10], -1) == 0
    assert bisect_right([0, 5, 10], 0) == 1
    assert bisect_right([0, 0, 0, 5, 10], 0) == 3
    assert bisect_right([0, 5, 10], 3) == 1
    assert bisect_right([0, 5, 10], 5) == 2
    assert bisect_right([0, 5, 5, 5, 10], 5) == 4
    assert bisect_right([0, 5, 10], 10) == 3
    assert bisect_right([0, 5, 10, 10, 10], 10) == 5
    assert bisect_right([0, 5, 10], 11) == 3


@pytest.mark.parametrize(
    "expected_tss_insertions, position, tss_vec, strand_vec",
    [
        ([0, 0, 1, 1, 1, 1, 0], 100, [5, 98, 99, 100, 101, 500], [1, 1, 1, 1, 1, -1]),
        ([0, 0, 0, 1, 2, 1, 0], 100, [5, 98, 99, 100, 101, 500], [1, 1, 1, 1, -1, -1]),
        ([0, 0, 0, 1, 3, 1, 0], 100, [5, 98, 99, 100, 101, 101, 500], [1, 1, 1, 1, -1, -1, -1]),
        ([0, 0, 1, 1, 2, 1, 0], 100, [5, 98, 99, 100, 101, 101, 500], [1, 1, 1, 1, -1, 1, -1]),
        ([0, 0, 0, 0, 0, 0, 0], -100, [5, 98, 99, 100, 101, 101, 500], [1, 1, 1, 1, -1, 1, -1]),
        ([0, 0, 0, 0, 0, 0, 0], 999, [5, 98, 99, 100, 101, 101, 500], [1, 1, 1, 1, -1, 1, -1]),
    ],
)
def test_update_tss_inserstions(
    expected_tss_insertions: list[int], position: int, tss_vec: list[int], strand_vec: list[int]
):
    """Verify that the update_tss_insertions function works as expected.

    - Use hand-calculated test cases.
    - Start from a random base to show that the tss_insertions are updated independently of their
    initial state.
    - Note that this code is optimized for speed and downstream of code that organizes tss_vec and
    strand_vec, so they are not tested for correctness. We only check with in-spec inputs here.
    """
    start_tss_insertions = np.random.randint(
        0,
        100,
        size=(
            len(
                expected_tss_insertions,
            )
        ),
    )
    actual_tss_insertions = update_tss_insertions(
        tss_insertions=start_tss_insertions,
        position=position,
        tss_vec=tss_vec,
        strand_vec=strand_vec,
    )
    assert (actual_tss_insertions - start_tss_insertions).tolist() == expected_tss_insertions
