import cython

DEF BOUNDS_CHECK = False


@cython.nogil
@cython.ccall
def bisect_left(arr: cython.int[:], v: cython.int) -> cython.size_t:
    """Cython thread-safe bisect left algorithm that operates on memory views."""
    with cython.boundscheck(BOUNDS_CHECK):
        delta: cython.size_t = len(arr)
        if delta == 0:
            return 0
        low: cython.size_t = 0
        while delta > 1:
            remainder: cython.size_t = delta % 2
            delta >>= 1
            check = low + delta
            if v > arr[check]:
                low = check
                delta += remainder
        return low if v <= arr[low] else low + delta

@cython.nogil
@cython.ccall
def bisect_right(arr: cython.int[:], v: cython.int) -> cython.size_t:
    """Cython thread-safe bisect right algorithm that operates on memory views."""
    with cython.boundscheck(BOUNDS_CHECK):
        delta: cython.size_t = len(arr)
        if delta == 0:
            return 0
        low: cython.size_t = 0
        while delta > 1:
            remainder: cython.size_t = delta % 2
            delta >>= 1
            check: cython.size_t = low + delta
            if v >= arr[check]:
                low = check
                delta += remainder
        return low if v < arr[low] else low + delta

@cython.nogil
@cython.ccall
@cython.inline
def update_insertions_range(
    tss_insertions: cython.ushort[:],
    tss_vec: cython.int[:],
    strand_vec: cython.int[:],
    start: cython.long,
    end: cython.long,
    half_window: cython.int,
    position: cython.int,
) -> cython.int:
    """Cython loop on memory views to update tss_insertions."""
    with cython.boundscheck(BOUNDS_CHECK):
        idx: cython.long
        for idx in range(start, end):
            tss: cython.int = tss_vec[idx]
            strand: cython.int = strand_vec[idx]
            tss_insertions[half_window + (position - tss) * strand] += 1
    return 0

@cython.nogil
@cython.ccall
def update_tss_insertions(
    tss_insertions: cython.ushort[:],
    position: cython.int,
    tss_vec: cython.int[:],
    strand_vec: cython.int[:],
) -> cython.int:
    """Update TSS counts with overlapping Transcription Start Site (TSS) in the given TSS vectors.

    Args:
        tss_insertions: Numpy array of counts of TSS.
        position: Genomic position to check (0-based).
        tss_vec: Numpy array of TSS positions on the same chromosome (0-based).
        strand_vec: Numpy array of strand signs corresponding to the TSS positions (+1 or -1).
    Returns:
        Numpy array of distances from each TSS to position (in the strand direction). If there are
        no TSS within half_window of the position, returns an empty array.
    """
    # find index where overlaps start on the left side
    half_window: cython.int = len(tss_insertions) // 2
    tss_left_idx: cython.long = bisect_left(tss_vec, position - half_window)
    # find relative (to left) index of end of overlaps on the right side
    tss_right_idx: cython.long = bisect_right(tss_vec, position + half_window)
    # update the insertions from the discovered range
    update_insertions_range(
        tss_insertions, tss_vec, strand_vec, tss_left_idx, tss_right_idx, half_window, position
    )
    return 0
