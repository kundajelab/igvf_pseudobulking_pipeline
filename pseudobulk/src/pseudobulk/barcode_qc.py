import dataclasses
from threading import Lock
from typing import Final

import numpy as np

from pseudobulk.fragment import Fragment
from pseudobulk.tss import update_tss_insertions
from pseudobulk.types import (
    COUNTS_ARRAY,
    COUNTS_DTYPE,
    POS_ARRAY,
    Barcode,
    Contig,
    PseudobulkName,
)

_MITO_CONTIG: Final[Contig] = Contig("chrM")
"""
Contig name for mitochondrial chromosome. Fragments on this chromosome will be ignored for TSS
enrichment calculations.
"""


@dataclasses.dataclass(slots=True, kw_only=True)
class BarcodeQc:
    """ATAC-seq QC statistics of one barcode, accumulated one fragment at a time."""

    tss_insertions: COUNTS_ARRAY
    fragments: list[Fragment] | None = dataclasses.field(default_factory=list)
    annotated: bool = False
    num_unique_frags: int = 0
    num_reads: int = 0
    mono_nucleosomal_frags: int = 0
    nucleosome_free_frags: int = 0
    lock: Lock = dataclasses.field(default_factory=Lock)

    @classmethod
    def new(cls, tss_half_window: int, is_pseudobulk: bool) -> BarcodeQc:
        """Create the QC of a barcode that has no fragments yet.

        Args:
            tss_half_window: Half the size of the window to check for TSS overlaps.
            is_pseudobulk: If True, the barcode is in a pseudobulk, so keep its fragments to be
                written out to the pseudobulk later.
        """
        return cls(
            tss_insertions=np.zeros((2 * tss_half_window + 1,), dtype=COUNTS_DTYPE),
            fragments=[] if is_pseudobulk else None,
        )

    def update_from_fragment(
        self,
        fragment: Fragment,
        tss_locs: dict[Contig, tuple[POS_ARRAY, POS_ARRAY]],
    ) -> None:
        """Update the QC with one more fragment of the barcode.

        Args:
            fragment: The fragment.
            tss_locs: Transcription Start Site (TSS) positions and strands on each contig, as loaded
                by utils.load_tss_locs.
        """
        shifted = fragment.shifted
        self.num_unique_frags += 1
        self.num_reads += fragment.num_reads
        frag_length = shifted.end - shifted.start
        self.nucleosome_free_frags += int(frag_length < 148)
        self.mono_nucleosomal_frags += int(148 <= frag_length < 295)
        if self.fragments is not None:
            self.fragments.append(fragment)  # this is a pseudobulk QC, store the fragment
        # TSS enrichment
        if (fragment.contig in tss_locs) and (fragment.contig != _MITO_CONTIG):
            tss_vec, strand_vec = tss_locs[fragment.contig]
            update_tss_insertions(
                self.tss_insertions,
                position=shifted.start,
                tss_vec=tss_vec,
                strand_vec=strand_vec,
            )
            update_tss_insertions(
                self.tss_insertions,
                position=shifted.end - 1,
                tss_vec=tss_vec,
                strand_vec=strand_vec,
            )

    @property
    def num_dup_reads(self) -> int:
        """Number of duplicate reads: the reads of each fragment after the first."""
        return self.num_reads - self.num_unique_frags

    @property
    def pct_duplicated_reads(self) -> float:
        """Percentage of reads that are duplicates."""
        return (self.num_dup_reads / self.num_reads) * 100

    @property
    def nucleosomal_signal(self) -> float:
        """Ratio of mono-nucleosomal to nucleosome-free fragments, adding 1 to each count."""
        return (1 + self.mono_nucleosomal_frags) / (1 + self.nucleosome_free_frags)

    @property
    def _tss_insertions_flank_mean(self) -> float:
        """Mean of TSS insertions at the flanks of the peak."""
        return (self.tss_insertions[:100].sum() + self.tss_insertions[-100:].sum()) / 200.0

    def _tss_insertions_center(self, tss_half_window: int, tss_half_smooth_window: int) -> float:
        """Mean of TSS insertions near the center of the peak."""
        start: int = tss_half_window - tss_half_smooth_window
        end: int = tss_half_window + tss_half_smooth_window + 1
        return self.tss_insertions[start:end].mean()

    def tss_enrichment(self, tss_half_window: int, tss_half_smooth_window: int) -> float:
        """Mean of TSS insertions near the center of the peak, relative to the flanks.

        Args:
            tss_half_window: Half the size of the window to check for TSS overlaps.
            tss_half_smooth_window: Half the size of the window to smooth over the peak of TSS
                insertions.
        """
        return (
            self._tss_insertions_center(
                tss_half_window=tss_half_window, tss_half_smooth_window=tss_half_smooth_window
            )
            / (
                self._tss_insertions_flank_mean + 0.1
            )  # add 0.1 like snapatac2 to avoid division by zero
        )

    @classmethod
    def header_columns(cls) -> list[str]:
        """Get the column names of the QC report, in the order of csv_columns."""
        return [
            "analysis_set_accession",
            "barcode_sample",
            "pseudobulk_id",
            "annotated",
            "num_frags",
            "pct_duplicated_reads",
            "nucleosomal_signal",
            "tss_enrichment",
            "raw-num_reads",
            "raw-num_dup_reads",
            "raw-mono_nucleosomal_frags",
            "raw-nucleosome_free_frags",
        ]

    def csv_columns(
        self,
        analysis_set_accession: str,
        barcode_sample: Barcode,
        pseudobulk_id: PseudobulkName | None,
        tss_half_window: int,
        tss_half_smooth_window: int,
    ) -> list[int | float | str]:
        """Get the values of this barcode's row of the QC report, in the order of header_columns.

        Args:
            analysis_set_accession: Accession of the analysis set the barcode is from.
            barcode_sample: The barcode.
            pseudobulk_id: ID of the pseudobulk the barcode is in, or None if it is in none.
            tss_half_window: Half the size of the window to check for TSS overlaps.
            tss_half_smooth_window: Half the size of the window to smooth over the peak of TSS
                insertions.
        """
        return [
            analysis_set_accession,
            barcode_sample,
            "null" if pseudobulk_id is None else pseudobulk_id,
            self.annotated,
            self.num_unique_frags,
            self.pct_duplicated_reads,
            self.nucleosomal_signal,
            self.tss_enrichment(
                tss_half_window=tss_half_window, tss_half_smooth_window=tss_half_smooth_window
            ),
            self.num_reads,
            self.num_dup_reads,
            self.mono_nucleosomal_frags,
            self.nucleosome_free_frags,
        ]
