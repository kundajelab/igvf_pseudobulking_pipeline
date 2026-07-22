import numpy as np

from pseudobulk.types import (
    COUNTS_ARRAY,
    COUNTS_DTYPE,
    POS_ARRAY,
    POS_DTYPE,
)

def bisect_left(arr: POS_ARRAY, pos: POS_DTYPE.type) -> np.int64: ...
def bisect_right(arr: POS_ARRAY, pos: POS_DTYPE.type) -> np.int64: ...
def update_insertions_range(
    tss_insertions: COUNTS_ARRAY,
    tss_vec: POS_ARRAY,
    strand_vec: POS_ARRAY,
    start: np.int64,
    end: np.int64,
    half_window: POS_DTYPE.type,
    position: POS_DTYPE.type,
) -> POS_DTYPE.type: ...
def update_tss_insertions(
    tss_insertions: COUNTS_DTYPE.type,
    position: POS_DTYPE.type,
    tss_vec: POS_ARRAY,
    strand_vec: POS_ARRAY,
) -> POS_DTYPE.type: ...
