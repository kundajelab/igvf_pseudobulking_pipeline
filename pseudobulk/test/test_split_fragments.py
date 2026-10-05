import gzip
from pathlib import Path
from threading import Thread
from typing import NamedTuple

import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from pseudobulk.tools.split_fragments import split_fragments

_ACCESSION = "ACC1"
_TSS_HALF_WINDOW = 2000

# (contig, start, end, barcode, num_reads)
_FRAGMENTS: list[tuple[str, int, int, str, int]] = [
    # shifted to 4904-5096 (192 bp, mono-nucleosomal), cutting either side of the chr1:5000 + TSS
    ("chr1", 4900, 5100, "bcA", 2),
    # shifted to 30004-30096 (92 bp, nucleosome-free), far from any TSS
    ("chr1", 30000, 30100, "bcA", 1),
    # on a contig missing from the chrom sizes: counted in QC, but not written to the pseudobulk
    ("chrUn_x", 100, 200, "bcA", 1),
    # shifted to 19954-20096 (142 bp, nucleosome-free), cutting either side of the chr1:20000 - TSS
    ("chr1", 19950, 20100, "bcB", 3),
    # shifted to 1004-1496 (492 bp, neither class)
    ("chr2", 1000, 1500, "bcC", 1),
    # a barcode that is not in the metadata
    ("chr1", 5000, 5050, "bcD", 1),
]


def _metadata_row(barcode: str, cell_name: str, accession: str) -> dict[str, str]:
    return {
        "barcode_sample": barcode,
        "cell_name": cell_name,
        "cell_description": f"a {cell_name}",
        "CL_id": f"CL:{cell_name}",
        "CL_term_name": cell_name,
        "subsample": "s0",
        "analysis_set_accession": accession,
    }


class _Inputs(NamedTuple):
    fragments_file: Path
    metadata_loc: Path
    chrom_sizes: Path
    tss_tsv: Path

    def split(self, output_dir: Path, num_threads: int) -> None:
        split_fragments(
            fragments_file=self.fragments_file,
            output_dir=output_dir,
            metadata_loc=self.metadata_loc,
            chrom_sizes=self.chrom_sizes,
            tss_tsv=self.tss_tsv,
            num_threads=num_threads,
        )


def _write_inputs(
    tmp_path: Path,
    fragments: list[tuple[str, int, int, str, int]],
    metadata_rows: list[dict[str, str]],
) -> _Inputs:
    fragments_file = tmp_path / f"{_ACCESSION}.bed"
    with open(fragments_file, "w") as f_out:
        for contig, start, end, barcode, num_reads in fragments:
            f_out.write(f"{contig}\t{start}\t{end}\t{barcode}\t{num_reads}\n")
    metadata = tmp_path / "metadata.tsv"
    pd.DataFrame(metadata_rows).to_csv(metadata, sep="\t", index=False)
    chrom_sizes = tmp_path / "chr_sizes.tsv"
    chrom_sizes.write_text("chr1\t100000\nchr2\t100000\n")
    tss_tsv = tmp_path / "tss.tsv"
    # TSS are 1-based in the TSV
    pd.DataFrame(
        {
            "gene": ["G1", "G2"],
            "transcript": ["G1-201", "G2-201"],
            "chro": ["chr1", "chr1"],
            "TSS": [5000, 20000],
            "strand": ["+", "-"],
        }
    ).to_csv(tss_tsv, sep="\t", index=False)
    return _Inputs(
        fragments_file=fragments_file,
        metadata_loc=metadata,
        chrom_sizes=chrom_sizes,
        tss_tsv=tss_tsv,
    )


def _read_lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


def _read_qc(output_dir: Path) -> tuple[pd.DataFrame, np.ndarray]:
    """Read the QC report, and the dense TSS matrix with rows in the same order."""
    qc_dir = output_dir / "atac_qc_reports"
    qc = pd.read_csv(qc_dir / f"{_ACCESSION}.tsv", sep="\t")
    tss = scipy.sparse.load_npz(qc_dir / f"{_ACCESSION}_tss_matrix.npz").toarray()
    return qc, tss


def test_split_fragments(tmp_path: Path) -> None:
    metadata_rows = [
        _metadata_row("bcA", "T cell", _ACCESSION),
        _metadata_row("bcB", "T cell", _ACCESSION),
        _metadata_row("bcC", "B cell", _ACCESSION),
        # annotated, but has no fragments
        _metadata_row("bcE", "T cell", _ACCESSION),
        # the same barcode in another accession must not affect this one
        _metadata_row("bcD", "NK cell", "ACC2"),
    ]
    inputs = _write_inputs(tmp_path, _FRAGMENTS, metadata_rows)
    output_dir = tmp_path / "output"

    inputs.split(output_dir=output_dir, num_threads=2)

    # pseudobulked fragments are the original fragments, restricted to the allowed contigs
    separated = output_dir / "separated_fragments"
    assert sorted(path.name for path in separated.iterdir()) == [
        f"B_cell-s0.{_ACCESSION}.tsv",
        f"T_cell-s0.{_ACCESSION}.tsv",
    ]
    assert sorted(_read_lines(separated / f"T_cell-s0.{_ACCESSION}.tsv")) == [
        "chr1\t19950\t20100\tbcB\t3",
        "chr1\t30000\t30100\tbcA\t1",
        "chr1\t4900\t5100\tbcA\t2",
    ]
    assert _read_lines(separated / f"B_cell-s0.{_ACCESSION}.tsv") == ["chr2\t1000\t1500\tbcC\t1"]

    # the pseudoreps hold the Tn5-shifted insertion sites, each one in exactly one of reps 1 and 2
    for pseudobulk, expected_insertions in [
        ("B_cell-s0", [("chr2", 1004), ("chr2", 1495)]),
        (
            "T_cell-s0",
            [
                ("chr1", 4904),
                ("chr1", 5095),
                ("chr1", 30004),
                ("chr1", 30095),
                ("chr1", 19954),
                ("chr1", 20095),
            ],
        ),
    ]:
        rep_t, rep_1, rep_2 = (
            _read_lines(
                output_dir / f"separated_pseudorep{rep}" / f"{pseudobulk}.{_ACCESSION}.{ext}.tsv"
            )
            for rep, ext in [("T", "t"), ("1", "1"), ("2", "2")]
        )
        insertions = sorted(
            (contig, int(start)) for contig, start, end, *_ in (line.split("\t") for line in rep_t)
        )
        assert insertions == sorted(expected_insertions)
        assert all(int(line.split("\t")[2]) == int(line.split("\t")[1]) + 1 for line in rep_t)
        assert sorted(rep_1 + rep_2) == sorted(rep_t)

    qc, tss = _read_qc(output_dir)
    qc = qc.set_index("barcode_sample")
    tss_rows = {barcode: tss[idx] for idx, barcode in enumerate(qc.index)}
    assert sorted(qc.index) == ["bcA", "bcB", "bcC", "bcD"]
    assert (qc["analysis_set_accession"] == _ACCESSION).all()
    assert qc["pseudobulk_id"].fillna("null").to_dict() == {
        "bcA": "T_cell-s0",
        "bcB": "T_cell-s0",
        "bcC": "B_cell-s0",
        "bcD": "null",
    }
    assert qc["annotated"].to_dict() == {"bcA": True, "bcB": True, "bcC": True, "bcD": False}
    assert qc["num_frags"].to_dict() == {"bcA": 3, "bcB": 1, "bcC": 1, "bcD": 1}
    assert qc["raw-num_reads"].to_dict() == {"bcA": 4, "bcB": 3, "bcC": 1, "bcD": 1}
    assert qc["raw-num_dup_reads"].to_dict() == {"bcA": 1, "bcB": 2, "bcC": 0, "bcD": 0}
    np.testing.assert_allclose(
        qc.loc[["bcA", "bcB", "bcC"], "pct_duplicated_reads"], [25, 200 / 3, 0]
    )
    assert qc["raw-mono_nucleosomal_frags"].to_dict() == {"bcA": 1, "bcB": 0, "bcC": 0, "bcD": 0}
    assert qc["raw-nucleosome_free_frags"].to_dict() == {"bcA": 2, "bcB": 1, "bcC": 0, "bcD": 1}
    np.testing.assert_allclose(
        qc.loc[["bcA", "bcB", "bcC"], "nucleosomal_signal"], [2 / 3, 1 / 2, 1]
    )

    # TSS insertion offsets are relative to the TSS, in the direction of its strand
    assert tss.shape == (4, 2 * _TSS_HALF_WINDOW + 1)
    assert np.flatnonzero(tss_rows["bcA"]).tolist() == [2000 - 95, 2000 + 96]
    assert np.flatnonzero(tss_rows["bcB"]).tolist() == [2000 - 96, 2000 + 45]
    assert not tss_rows["bcC"].any()
    # bcD is shifted to 5004-5046, so cuts at 5004 and 5045, 5 and 46 bp downstream of the TSS
    assert np.flatnonzero(tss_rows["bcD"]).tolist() == [2005, 2046]
    # the center is the 11 bins around the TSS, and 0.1 is added to the (empty) flanks
    np.testing.assert_allclose(qc.loc[["bcA", "bcB", "bcC"], "tss_enrichment"], 0.0)
    np.testing.assert_allclose(qc.loc["bcD", "tss_enrichment"], (1 / 11) / 0.1)


def _random_fragments(
    num_barcodes: int, num_per_barcode: int
) -> list[tuple[str, int, int, str, int]]:
    rng = np.random.default_rng(seed=0)
    fragments: list[tuple[str, int, int, str, int]] = []
    for idx in range(num_barcodes * num_per_barcode):
        start = int(rng.integers(0, 90_000))
        fragments.append(
            (
                "chr1",
                start,
                start + int(rng.integers(20, 600)),
                f"bc{idx % num_barcodes}",
                int(rng.integers(1, 4)),
            )
        )
    # fragments files are sorted by position, which interleaves the barcodes
    return sorted(fragments)


def test_split_fragments_is_reproducible(tmp_path: Path) -> None:
    """Outputs must be identical whatever the number of threads, and so the order they ran in."""
    num_barcodes = 20
    metadata_rows = [
        _metadata_row(f"bc{idx}", "T cell" if idx % 2 else "B cell", _ACCESSION)
        for idx in range(num_barcodes)
    ]
    inputs = _write_inputs(tmp_path, _random_fragments(num_barcodes, 2_500), metadata_rows)

    output_dirs: list[Path] = []
    for num_threads in (1, 4):
        output_dir = tmp_path / f"output_{num_threads}"
        inputs.split(output_dir=output_dir, num_threads=num_threads)
        output_dirs.append(output_dir)

    def output_files(output_dir: Path) -> dict[Path, bytes]:
        return {
            path.relative_to(output_dir): path.read_bytes()
            for path in sorted(output_dir.rglob("*.tsv"))
        }

    single_thread, multi_thread = (output_files(output_dir) for output_dir in output_dirs)
    assert list(single_thread) == list(multi_thread)
    assert len(single_thread) == 4 * 2 + 1
    for path, contents in single_thread.items():
        assert contents == multi_thread[path], f"{path} differs between thread counts"
    # the TSS matrix is a zip file holding modification times, so compare its contents instead
    single_tss, multi_tss = (_read_qc(output_dir)[1] for output_dir in output_dirs)
    np.testing.assert_array_equal(single_tss, multi_tss)


def test_split_fragments_fails_on_malformed_line(tmp_path: Path) -> None:
    """A worker thread that fails must fail the whole run, rather than only losing its fragments."""
    metadata_rows = [_metadata_row("bcA", "T cell", _ACCESSION)]
    inputs = _write_inputs(tmp_path, _FRAGMENTS[:2], metadata_rows)
    with open(inputs.fragments_file, "a") as f_out:
        f_out.write("chr1\t100\n")  # truncated line
    with pytest.raises(ValueError, match="not enough values to unpack"):
        inputs.split(output_dir=tmp_path / "output", num_threads=1)


def test_split_fragments_fails_on_truncated_gzip(tmp_path: Path) -> None:
    """The reader thread failing must fail the run, rather than look like the end of the file."""
    metadata_rows = [_metadata_row(f"bc{idx}", "T cell", _ACCESSION) for idx in range(4)]
    inputs = _write_inputs(tmp_path, _random_fragments(4, 5_000), metadata_rows)
    fragments_gz = tmp_path / f"{_ACCESSION}.bed.gz"
    compressed = gzip.compress(inputs.fragments_file.read_bytes())
    fragments_gz.write_bytes(compressed[: len(compressed) // 2])
    with pytest.raises(EOFError):
        inputs._replace(fragments_file=fragments_gz).split(
            output_dir=tmp_path / "output", num_threads=2
        )


def test_split_fragments_fails_when_every_worker_fails(tmp_path: Path) -> None:
    """With every worker thread failed, the reader must stop too.

    Otherwise it would wait forever for the workers to make room in the queue.
    """
    num_threads = 2
    metadata_rows = [_metadata_row(f"bc{idx}", "T cell", _ACCESSION) for idx in range(4)]
    # more fragments than the queue holds, so that the reader has to wait for the workers
    inputs = _write_inputs(tmp_path, _random_fragments(4, 10_000), metadata_rows)
    lines = inputs.fragments_file.read_text()
    inputs.fragments_file.write_text("bad line\n" * num_threads + lines)

    # run in a daemon thread, so that a hang fails the test rather than blocking the test run
    exceptions: list[BaseException] = []

    def _split() -> None:
        try:
            inputs.split(output_dir=tmp_path / "output", num_threads=num_threads)
        except BaseException as exception:
            exceptions.append(exception)

    thread = Thread(target=_split, daemon=True)
    thread.start()
    thread.join(timeout=60)
    assert not thread.is_alive(), "split_fragments hung after its worker threads failed"
    assert len(exceptions) == 1
    assert isinstance(exceptions[0], ValueError)
