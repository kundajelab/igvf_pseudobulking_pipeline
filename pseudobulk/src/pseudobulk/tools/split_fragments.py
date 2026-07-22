import csv
import dataclasses
import gzip
import itertools
import logging
import random
import time
from collections import deque
from collections.abc import Mapping
from contextlib import ExitStack
from io import BufferedReader
from pathlib import Path
from threading import (
    Lock,
    Thread,
    local,
)
from typing import (
    BinaryIO,
    TextIO,
)

import numpy as np
import pandas as pd
import psutil
import scipy.sparse

from pseudobulk import utils
from pseudobulk.barcode_qc import BarcodeQc
from pseudobulk.fragment import Fragment
from pseudobulk.types import (
    COUNTS_DTYPE,
    POS_ARRAY,
    Barcode,
    Contig,
    PseudobulkName,
)
from pseudobulk.utils import create_and_write


@dataclasses.dataclass
class PseudobulkFiles:
    pseudorep1_out: TextIO
    pseudorep2_out: TextIO
    pseudorep_t_out: TextIO
    fragments_out: TextIO

    @classmethod
    def new(
        cls,
        pseudobulk: PseudobulkName,
        exit_stack: ExitStack,
        output_dir: Path,
        analysis_set_accession: str,
    ) -> "PseudobulkFiles":
        file_base = f"{pseudobulk}"
        return PseudobulkFiles(
            pseudorep1_out=exit_stack.enter_context(
                create_and_write(
                    output_dir
                    / "separated_pseudorep1"
                    / f"{file_base}.{analysis_set_accession}.1.tsv",
                    mode="at",
                )
            ),
            pseudorep2_out=exit_stack.enter_context(
                create_and_write(
                    output_dir
                    / "separated_pseudorep2"
                    / f"{file_base}.{analysis_set_accession}.2.tsv",
                    mode="at",
                )
            ),
            pseudorep_t_out=exit_stack.enter_context(
                create_and_write(
                    output_dir
                    / "separated_pseudorepT"
                    / f"{file_base}.{analysis_set_accession}.t.tsv",
                    mode="at",
                )
            ),
            fragments_out=exit_stack.enter_context(
                create_and_write(
                    output_dir
                    / "separated_fragments"
                    / f"{file_base}.{analysis_set_accession}.tsv",
                    mode="at",
                )
            ),
        )

    def write_fragment(self, fragment: Fragment) -> None:
        shifted = fragment.shifted
        start_point = shifted.start_point
        end_point = shifted.end_point
        self.pseudorep_t_out.write(f"{start_point}")
        self.pseudorep_t_out.write(f"{end_point}")
        if random.random() < 0.5:
            self.pseudorep1_out.write(f"{start_point}")
            self.pseudorep1_out.write(f"{end_point}")
        else:
            self.pseudorep2_out.write(f"{start_point}")
            self.pseudorep2_out.write(f"{end_point}")
        self.fragments_out.write(f"{fragment}")


@dataclasses.dataclass(slots=True, kw_only=True)
class SharedThreadData:
    fragments_in: BinaryIO
    fragments_deque: deque[bytes | None]
    fragments_in_lock: Lock
    thread_local_data: local
    tss_half_window: int
    tss_locs: dict[Contig, tuple[POS_ARRAY, POS_ARRAY]]
    pseudobulk_barcodes: set[Barcode]
    barcode_qcs: dict[Barcode, BarcodeQc]
    barcode_qcs_lock: Lock
    output_dir: Path
    shutdown_lock: Lock
    logger: logging.Logger
    num_workers: int
    num_lines: int = 0
    num_shutdown: int = 0

    @classmethod
    def new(
        cls,
        fragments_in: BinaryIO,
        fragments_deque: deque[bytes | None],
        tss_half_window: int,
        tss_tsv: Path,
        pseudobulk_barcodes: set[Barcode],
        output_dir: Path,
        logger: logging.Logger,
        num_workers: int,
    ) -> "SharedThreadData":
        return cls(
            fragments_in=fragments_in,
            fragments_deque=fragments_deque,
            fragments_in_lock=Lock(),
            thread_local_data=local(),
            tss_half_window=tss_half_window,
            tss_locs=utils.load_tss_locs(tss_tsv),
            pseudobulk_barcodes=pseudobulk_barcodes,
            barcode_qcs={},
            barcode_qcs_lock=Lock(),
            output_dir=output_dir,
            shutdown_lock=Lock(),
            logger=logger,
            num_workers=num_workers,
        )

    def _barcode_qc_factory(self, barcode_sample: Barcode) -> BarcodeQc:
        return BarcodeQc.new(
            tss_half_window=self.tss_half_window,
            is_pseudobulk=barcode_sample in self.pseudobulk_barcodes,
        )

    @property
    def index_buffer(self) -> POS_ARRAY:
        index_buffer = getattr(self.thread_local_data, "index_buffer", None)
        if index_buffer is None:
            index_buffer = np.empty((2 * self.tss_half_window + 1,), dtype=utils.POS_DTYPE)
            self.thread_local_data.index_buffer = index_buffer
        return index_buffer

    @property
    def next_fragment(self) -> Fragment | None:
        while True:
            try:
                if len(self.fragments_deque) > 0:
                    fragment_line = self.fragments_deque.popleft()
                    break
                else:
                    time.sleep(1e-6)
            except IndexError:
                time.sleep(1e-6)
        if fragment_line is None:
            return None
        self.num_lines += 1
        return Fragment.from_line(fragment_line.decode("utf-8"))

    def get_barcode_qc(self, barcode_sample: Barcode) -> BarcodeQc:
        barcode_qc = self.barcode_qcs.get(barcode_sample, None)
        if barcode_qc is None:
            with self.barcode_qcs_lock:
                # double-check that a different thread didn't already create this BarcodeQc while we
                # waited to acquire the lock
                barcode_qc = self.barcode_qcs.get(barcode_sample, None)
                if barcode_qc is None:
                    # this Barcode is new, create a new BarcodeQc
                    barcode_qc = self._barcode_qc_factory(barcode_sample)
                    self.barcode_qcs[barcode_sample] = barcode_qc
        return barcode_qc

    def process_fragment(self, fragment: Fragment) -> None:
        # get the correct BarcodeQc
        barcode_qc = self.get_barcode_qc(fragment.barcode_sample)

        # update that BarcodeQc with data from the fragment
        with barcode_qc.lock:
            barcode_qc.update_from_fragment(
                fragment=fragment,
                tss_locs=self.tss_locs,
            )

    def log_progress(self, elapsed_time: float) -> None:
        human_elapsed_time = utils.elapsed_time(elapsed_time)
        self.logger.info(
            f"Processed {self.num_lines} lines in {human_elapsed_time} "
            f"({self.num_lines / elapsed_time:.1f} lines/s)"
        )
        self.logger.info(f"Have {len(self.barcode_qcs)} unique barcode QCs")


def _thread_qc(shared_thread_data: SharedThreadData) -> None:
    try:
        fragment = shared_thread_data.next_fragment
        while fragment is not None:
            shared_thread_data.process_fragment(fragment)
            fragment = shared_thread_data.next_fragment
    finally:
        with shared_thread_data.shutdown_lock:
            shared_thread_data.num_shutdown += 1


def _fill_deque(shared_thread_data: SharedThreadData) -> None:
    try:
        batch_size = 5000
        max_queue_size = 10000
        sleep_time = 0.001
        fragments_deque = shared_thread_data.fragments_deque

        for batch in itertools.batched(shared_thread_data.fragments_in, batch_size):
            while len(fragments_deque) > max_queue_size:
                time.sleep(sleep_time)
            fragments_deque.extend(batch)

    finally:
        for _ in range(shared_thread_data.num_workers):
            shared_thread_data.fragments_deque.append(None)
        with shared_thread_data.shutdown_lock:
            shared_thread_data.num_shutdown += 1


def _run_qc(
    fragments_file: Path,
    num_workers: int,
    tss_half_window: int,
    tss_tsv: Path,
    output_dir: Path,
    pseudobulk_barcodes: set[Barcode],
    logger: logging.Logger,
) -> dict[Barcode, BarcodeQc]:
    opener = gzip.open if fragments_file.suffix == ".gz" else open
    with (
        opener(f"{fragments_file}", "rb") as fragments_in_raw,
        BufferedReader(fragments_in_raw) as fragments_in,
    ):
        fragments_deque: deque[bytes | None] = deque()
        shared_thread_data = SharedThreadData.new(
            fragments_in=fragments_in,
            fragments_deque=fragments_deque,
            tss_half_window=tss_half_window,
            tss_tsv=tss_tsv,
            pseudobulk_barcodes=pseudobulk_barcodes,
            output_dir=output_dir,
            logger=logger,
            num_workers=num_workers,
        )
        threads = tuple(
            Thread(target=_thread_qc, name=f"Thread-{idx}", args=(shared_thread_data,))
            for idx in range(num_workers)
        ) + (Thread(target=_fill_deque, name=f"Thread-{num_workers}", args=(shared_thread_data,)),)
        for thread in threads:
            thread.start()
        start_time = time.time()
        last_time = start_time
        while shared_thread_data.num_shutdown < num_workers + 1:
            time.sleep(1.0)
            current_time = time.time()
            if current_time - last_time > 10:
                shared_thread_data.log_progress(elapsed_time=current_time - start_time)
                last_time = current_time
        logger.info(
            f"Finished processing {shared_thread_data.num_lines} lines, waiting for thread pool"
            " to shut down."
        )
        for thread in threads:
            thread.join()
    shared_thread_data.log_progress(elapsed_time=time.time() - start_time)
    return shared_thread_data.barcode_qcs


def _write_qc(
    analysis_set_accession: str,
    barcode_qcs: dict[Barcode, BarcodeQc],
    barcodes_to_pseudobulks: Mapping[Barcode, PseudobulkName],
    output_dir: Path,
    tss_half_window: int,
    tss_half_smooth_window: int,
    logger: logging.Logger,
) -> None:
    """Write QC reports and sparse matrix of Transcription Start Sites."""
    logger.info("Writing QC reports")
    barcodes = list(barcode_qcs.keys())
    num_barcodes = len(barcodes)
    tss_row_len = 2 * tss_half_window + 1
    # Fill the compressed row buffers of the TSS matrix directly, instead of building one sparse array
    # per barcode and stacking them at the end: scipy.sparse.vstack has to allocate a second copy of
    # every insertion while the per-barcode arrays are still alive, which doubles peak memory. Each
    # per-barcode sparse array also costs ~1 KB of container overhead on top of its insertions, which
    # dominates for the many barcodes that have very few insertions.
    # NOTE: the buffers are sized exactly, using one cheap pass to count the non-zero insertions, so
    # that they never have to be grown (which would copy them as well).
    total_insertions = sum(
        int(np.count_nonzero(barcode_qcs[barcode_sample].tss_insertions))
        for barcode_sample in barcodes
    )
    # NOTE: the indices and the row pointers must share a dtype, otherwise scipy copies the narrower
    # of the two to widen it when the matrix is created.
    index_dtype = np.int64 if total_insertions > 2**31 - 1 else np.int32
    tss_insertion_counts = np.empty((total_insertions,), dtype=COUNTS_DTYPE)
    tss_insertion_locs = np.empty((total_insertions,), dtype=index_dtype)
    tss_row_starts = np.zeros((num_barcodes + 1,), dtype=index_dtype)
    num_written = 0
    with utils.create_and_write(
        output_dir / "atac_qc_reports" / f"{analysis_set_accession}.tsv"
    ) as qc_out:
        writer = csv.writer(qc_out, delimiter="\t")
        writer.writerow(BarcodeQc.header_columns())
        for idx, barcode_sample in enumerate(barcodes):
            barcode_qc = barcode_qcs.pop(barcode_sample)
            writer.writerow(
                barcode_qc.csv_columns(
                    analysis_set_accession=analysis_set_accession,
                    barcode_sample=barcode_sample,
                    pseudobulk_id=barcodes_to_pseudobulks.get(barcode_sample, None),
                    tss_half_window=tss_half_window,
                    tss_half_smooth_window=tss_half_smooth_window,
                )
            )
            (insertion_locs,) = np.nonzero(barcode_qc.tss_insertions)
            num_insertions = insertion_locs.size
            insertions_slice = slice(num_written, num_written + num_insertions)
            tss_insertion_locs[insertions_slice] = insertion_locs
            tss_insertion_counts[insertions_slice] = barcode_qc.tss_insertions[insertion_locs]
            num_written += num_insertions
            tss_row_starts[idx + 1] = num_written
            del barcode_qc.tss_insertions
            del barcode_qc
            if idx % 10000 == 9999:
                barcode_qcs = {key: val for key, val in barcode_qcs.items()}

    del barcode_qcs
    logger.info("Writing TSS sparse matrix.")
    # NOTE: creating the matrix from the buffers does not copy them
    scipy.sparse.save_npz(
        output_dir / "atac_qc_reports" / f"{analysis_set_accession}_tss_matrix.npz",
        scipy.sparse.csr_array(
            (tss_insertion_counts, tss_insertion_locs, tss_row_starts),
            shape=(num_barcodes, tss_row_len),
        ),
    )


def split_fragments(
    *,
    fragments_file: Path,
    output_dir: Path,
    metadata_loc: Path,
    chrom_sizes: Path,
    tss_tsv: Path,
    tss_half_window: int = 2000,
    tss_half_smooth_window: int = 5,
    num_threads: int = 0,
    random_seed: int = 42,
):
    """
    Process a fragments file and split it into pseudobulks.

    Args:
        fragments_file: Path to fragments (.bed.gz) file to split.
        output_dir: Directory to write output files to.
        metadata_loc: Path to metadata file. The file should be a tab-separated values (TSV) file
            with columns "analysis_set_accession", "barcode_sample", "subsample", "subsample_orig", and
            "annotation"
        chrom_sizes: Path to TSV with contig sizes
        tss_tsv: Path to TSV with species-dependent Transcription Start Sites (TSS)
        tss_half_window: half_window: Half the size of the window to check for TSS overlaps.
        tss_half_smooth_window: Half the size of the window to smooth over the peak of TSS
            insertions.
        num_threads: Number of threads to use for processing. If <= 0, will use the number of
            physical CPUs available.
        random_seed: Arbitrary random number to generate deterministic results.
    """
    logger = logging.getLogger(name=f"{__package__} split-fragments")

    # Load metadata
    metadata_df = utils.load_metadata(metadata_loc)
    logger.info(f"Loaded metadata with {len(metadata_df)} rows")
    # Subset metadata to current analysis accession (fragment file name)
    analysis_set_accession = fragments_file.name.split(".")[0]
    metadata_df: pd.DataFrame = metadata_df[
        metadata_df["analysis_set_accession"] == f"{analysis_set_accession}"
    ]
    logger.info(
        f"Subset metadata to {len(metadata_df)} rows with for accession {analysis_set_accession}"
    )

    # Compute barcodes --> annotation mapping
    barcodes_to_pseudobulks = utils.map_barcodes_to_pseudobulks(metadata_df=metadata_df)

    # Get allowed chromosomes
    chrom_sizes_df = utils.read_csv(chrom_sizes, names=["chr", "size"])
    allowed_chrs: set[Contig] = set(chrom_sizes_df["chr"].unique().tolist())

    # Set the maximum number of CPUs to use for processing
    num_workers: int = psutil.cpu_count(logical=False) if num_threads <= 0 else num_threads
    logger.info(f"Processing fragments file {fragments_file} with {num_workers} workers...")

    # Iterate through fragments file, updating a BarcodeQc object for each Barcode
    barcode_qcs = _run_qc(
        fragments_file=fragments_file,
        num_workers=num_workers,
        tss_half_window=tss_half_window,
        tss_tsv=tss_tsv,
        output_dir=output_dir,
        pseudobulk_barcodes=set(barcodes_to_pseudobulks.keys()),
        logger=logger,
    )

    logger.info("Updating pseudobulk stats")
    random.seed(random_seed)
    num_pseudobulk_barcodes: int = 0

    # all the fragments are QC-ed and split, write fragments to appropriate pseudobulk files
    for barcode_sample, barcode_qc in barcode_qcs.items():
        if barcode_qc.fragments is None:
            continue  # this is not a pseuobulk QC

        pseudobulk: PseudobulkName | None = barcodes_to_pseudobulks.get(barcode_sample, None)
        if pseudobulk is None:
            continue

        num_pseudobulk_barcodes += 1

        with ExitStack() as exit_stack:
            pseudobulk_out = PseudobulkFiles.new(
                pseudobulk=pseudobulk,
                exit_stack=exit_stack,
                output_dir=output_dir,
                analysis_set_accession=analysis_set_accession,
            )
            for fragment in barcode_qc.fragments:
                # Skip nonstandard chromosomes
                if fragment.contig in allowed_chrs:
                    barcode_qc.annotated = True
                    pseudobulk_out.write_fragment(fragment)
            barcode_qc.fragments = []

    logger.info(f"Wrote {num_pseudobulk_barcodes} pseudobulk barcodes.")

    _write_qc(
        analysis_set_accession=analysis_set_accession,
        barcode_qcs=barcode_qcs,
        barcodes_to_pseudobulks=barcodes_to_pseudobulks,
        output_dir=output_dir,
        tss_half_window=tss_half_window,
        tss_half_smooth_window=tss_half_smooth_window,
        logger=logger,
    )
