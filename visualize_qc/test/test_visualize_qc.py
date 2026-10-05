import gzip
import logging
import math
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.offline
import polars as pl
import pytest

from visualize_qc.visualize_qc import (
    ACCESSION_QC_SCHEMA,
    PSEUDOBULK_QC_SCHEMA,
    DataScales,
    ScaleInfo,
    _add_figure,
    _find_qc_files,
    _fmt_text,
    _get_outliers,
    _strip_trailing_zeros,
    _visualize_categorical,
    _visualize_numeric,
    _visualize_table,
    scan_csv,
    visualize_qc,
)

LOGGER = logging.getLogger("test-visualize-qc")

RNA_COLS = ("rna_read_count", "gene_count", "pct_mito", "pct_ribo")
ATAC_COLS = (
    "num_frags",
    "pct_duplicated_reads",
    "nucleosomal_signal",
    "tss_enrichment",
)
FRONT_COLS = (
    "analysis_set_accession",
    "barcode_sample",
    "annotated",
    "found_in_rna",
    "found_in_atac",
    "pseudobulk_id",
)
# Column orders written by pseudobulk combine-accession-qc: ATAC-only accessions put the ATAC
# columns before the (empty) RNA columns.
ACCESSION_COLS_WITH_RNA = (*FRONT_COLS, *RNA_COLS, *ATAC_COLS)
ACCESSION_COLS_ATAC_ONLY = (*FRONT_COLS, *ATAC_COLS, *RNA_COLS)
PSEUDOBULK_COLS = tuple(PSEUDOBULK_QC_SCHEMA.names())
TABLE_COLS = (
    "accession",
    "frag",
    "matrix",
    "only-frag",
    "only-matrix",
    "both",
    "neither",
)


def _write_tsv(path: Path, columns: Sequence[str], rows: Sequence[Sequence[object]]) -> Path:
    """Write a TSV (gzipped if the suffix is .gz), writing None as empty like pandas does."""
    lines = ["\t".join(columns)]
    lines.extend("\t".join("" if val is None else f"{val}" for val in row) for row in rows)
    text = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as gz_file:
            gz_file.write(text)
    else:
        path.write_text(text)
    return path


def _accession_row(
    barcode: str, num_frags: int | None, rna_read_count: int | None
) -> dict[str, object]:
    return {
        "analysis_set_accession": "IGVFDS0000AAAA",
        "barcode_sample": barcode,
        "annotated": True,
        "found_in_rna": rna_read_count is not None,
        "found_in_atac": num_frags is not None,
        "pseudobulk_id": "pb1",
        "rna_read_count": rna_read_count,
        "gene_count": None if rna_read_count is None else rna_read_count // 2,
        "pct_mito": None if rna_read_count is None else 1.5,
        "pct_ribo": None if rna_read_count is None else 2.5,
        "num_frags": num_frags,
        "pct_duplicated_reads": None if num_frags is None else 0.25,
        "nucleosomal_signal": None if num_frags is None else 0.75,
        "tss_enrichment": None if num_frags is None else "inf",
    }


def _write_dict_rows(path: Path, columns: Sequence[str], rows: list[dict[str, object]]) -> Path:
    return _write_tsv(path, columns, [[row[col] for col in columns] for row in rows])


def test_scan_csv_matches_schema_by_name(tmp_path: Path) -> None:
    rows = [_accession_row("bc1", num_frags=10, rna_read_count=None)]
    path = _write_dict_rows(tmp_path / "a_per_cell_qc.tsv.gz", ACCESSION_COLS_ATAC_ONLY, rows)

    qc_df, index = scan_csv(path, index_col="barcode_sample", schema=ACCESSION_QC_SCHEMA)
    collected = qc_df.collect()

    expected_names = [col for col in ACCESSION_QC_SCHEMA.names() if col != "barcode_sample"]
    assert collected.columns == expected_names
    assert collected.schema == pl.Schema({col: ACCESSION_QC_SCHEMA[col] for col in expected_names})
    assert index.to_list() == ["bc1"]
    assert collected["num_frags"].to_list() == [10]
    assert collected["pct_duplicated_reads"].to_list() == [0.25]
    assert collected["rna_read_count"].to_list() == [None]
    assert collected["found_in_rna"].to_list() == [False]
    assert collected["tss_enrichment"].to_list() == [math.inf]


def test_scan_csv_rejects_mismatched_columns(tmp_path: Path) -> None:
    columns = ACCESSION_COLS_WITH_RNA[:-1]
    path = _write_tsv(tmp_path / "a.tsv", columns, [["x"] * len(columns)])
    with pytest.raises(ValueError, match="tss_enrichment"):
        scan_csv(path, schema=ACCESSION_QC_SCHEMA)


def test_scan_csv_does_not_mutate_exclude_cols(tmp_path: Path) -> None:
    path = _write_tsv(tmp_path / "a.tsv", ["idx", "a", "b"], [["x", 1, 2]])
    exclude_cols = {"b"}
    qc_df, _index = scan_csv(path, index_col="idx", exclude_cols=exclude_cols)
    assert exclude_cols == {"b"}
    assert qc_df.collect_schema().names() == ["a"]


def test_scan_csv_unknown_suffix(tmp_path: Path) -> None:
    path = _write_tsv(tmp_path / "a.txt", ["a"], [[1]])
    with pytest.raises(ValueError, match="Could not infer separator"):
        scan_csv(path)


def test_data_scales_from_data() -> None:
    data = pl.Series("x", [1.0, 2.0, 3.0, 4.0, 5.0, math.inf])
    scales = DataScales.from_data(data)
    assert scales.min_val == 1.0
    assert scales.max_val == math.inf
    assert scales.max_finite_val == 5.0
    assert scales.median == 3.0
    assert scales.iqr == scales.quartile_3 - scales.quartile_1


def test_data_scales_empty() -> None:
    scales = DataScales.from_data(pl.Series("x", [], dtype=pl.Float64))
    assert all(math.isnan(val) for val in scales.__dict__.values())
    assert ScaleInfo.from_data_scales(scales).unplotable_range


@pytest.mark.parametrize(
    ("has_minus_inf", "has_plus_inf", "expected_range"),
    [
        (False, False, (0.0, 10.0)),
        (True, False, (-5.0, 10.0)),
        (False, True, (0.0, 15.0)),
        (True, True, (-5.0, 15.0)),
    ],
)
def test_scale_info_range_extends_for_infinities(
    has_minus_inf: bool, has_plus_inf: bool, expected_range: tuple[float, float]
) -> None:
    scale_info = ScaleInfo(
        min_val=0.0,
        max_val=10.0,
        has_minus_inf=has_minus_inf,
        has_plus_inf=has_plus_inf,
    )
    assert scale_info.range == expected_range


@pytest.mark.parametrize("val", [-1234.5, -0.5, 0.0, 0.0005, 0.5, 1.0, 7.0, 1e6])
def test_symlog_inverts(val: float) -> None:
    scale_info = ScaleInfo(min_val=0.0, max_val=1.0, log_scale=0.5)
    assert scale_info.inv_scale(scale_info.scale(val)) == pytest.approx(val)


def test_symlog_ticks_are_consistent() -> None:
    data = pl.Series("x", [0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1000.0])
    scale_info = ScaleInfo.from_data_scales(DataScales.from_data(data))
    assert scale_info.log_scale is not None
    ticktext = scale_info.ticktext
    tickvals = scale_info.tickvals
    assert ticktext is not None
    assert tickvals is not None
    assert ticktext == ["0", "0.5", "1", "10", "100", "1000"]
    # positions are exact, not derived from rounded labels
    assert tickvals == [scale_info.scale(float(label)) for label in ticktext]


@pytest.mark.parametrize(
    ("min_val", "max_val", "log_scale"),
    [
        (0.0, 100.0, 1.0),
        (0.0, 1000.0, 1e-3),
        (0.002, 0.8, 1e-3),
        (-50.0, 2000.0, 0.5),
        (0.0, 1e15, 1e-3),
        (-1e20, 0.0, 1.0),
    ],
)
def test_symlog_ticks_are_round_and_not_crowded(
    min_val: float, max_val: float, log_scale: float
) -> None:
    scale_info = ScaleInfo(min_val=min_val, max_val=max_val, log_scale=log_scale)
    ticktext = scale_info.ticktext
    tickvals = scale_info.tickvals
    assert ticktext is not None
    assert tickvals is not None
    assert len(ticktext) == len(tickvals) >= 2
    values = [float(label) for label in ticktext]
    assert all(min_val <= val <= max_val for val in values)
    # every tick is 0 or 1, 2, or 5 times a power of 10
    assert all(
        val == 0 or f"{abs(val):e}".split("e")[0] in {"1.000000", "2.000000", "5.000000"}
        for val in values
    )
    bottom, top = scale_info.transformed_range
    min_gap = (top - bottom) / (2 * scale_info.num_ticks)
    assert np.all(np.diff(tickvals) >= min_gap)


@pytest.mark.parametrize(
    ("min_val", "max_val", "log_scale", "expected"),
    [
        (0.0, 100.0, 1.0, ["0", "0.5", "1", "5", "10", "20", "50", "100"]),
        (
            0.0,
            1e15,
            1e-3,
            ["0", "50", "1000", "100000", "10000000", "1e+09", "1e+11", "1e+13", "1e+15"],
        ),
        # no round numbers in range, so use the ends of the range
        (5.9, 9.0, 1.0, ["5.9", "9"]),
    ],
)
def test_symlog_tick_labels(
    min_val: float, max_val: float, log_scale: float, expected: list[str]
) -> None:
    scale_info = ScaleInfo(min_val=min_val, max_val=max_val, log_scale=log_scale)
    assert scale_info.ticktext == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("100", "100"),
        ("0.500", "0.5"),
        ("-5.900", "-5.9"),
        ("1.000e+20", "1e+20"),
        ("<1e-3", "<1e-3"),
    ],
)
def test_strip_trailing_zeros(text: str, expected: str) -> None:
    assert _strip_trailing_zeros(text) == expected


@pytest.mark.parametrize(
    ("val", "expected"),
    [
        (3.0, "3"),
        (0.125, "0.125"),
        (12345678.0, "12345678"),
        (123456789.0, "1.235e+08"),
        (1.5e20, "1.500e+20"),
        (0.0004, "<1e-3"),
        (-0.0004, ">-1e-3"),
        (0.0, "0"),
        (math.nan, "nan"),
        (math.inf, "inf"),
        (-math.inf, "-inf"),
    ],
)
def test_fmt_text(val: float, expected: str) -> None:
    assert _fmt_text(val) == expected


def test_get_outliers_samples_evenly() -> None:
    values = pl.Series("x", [-100, *range(10), 100, 200, 300, 400], dtype=pl.Int64)
    index = pl.Series("idx", [f"bc{i}" for i in range(len(values))])
    scale_info = ScaleInfo(min_val=-100, max_val=400)
    positions, outliers, outlier_index = _get_outliers(
        values, index, low=0, high=9, scale_info=scale_info, max_outliers=3
    )
    assert positions.to_list() == [-100, 200, 400]
    assert outliers.to_list() == [-100, 200, 400]
    assert outlier_index.to_list() == ["bc0", "bc12", "bc14"]


def test_infinite_outliers_drawn_at_range_edges() -> None:
    values = pl.Series("x", [-math.inf, *range(10), 100.0, math.inf], dtype=pl.Float64)
    index = pl.Series("idx", [f"bc{i}" for i in range(len(values))])
    scale_info, traces = _visualize_numeric(values, index, logger=LOGGER)
    assert scale_info.log_scale is not None
    bottom, top = scale_info.transformed_range
    assert bottom < scale_info.scale(0.0) < scale_info.scale(100.0) < top
    outlier_trace = list(traces)[2]
    assert isinstance(outlier_trace, go.Scattergl)
    assert outlier_trace.y is not None
    assert outlier_trace.hovertext is not None
    assert list(outlier_trace.y) == [bottom, scale_info.scale(100.0), top]
    assert list(outlier_trace.hovertext) == ["bc0", "bc11", "bc12"]
    assert outlier_trace.customdata is not None
    assert list(outlier_trace.customdata) == ["-inf", "100", "inf"]


def test_nan_is_dropped_not_an_outlier() -> None:
    values = pl.Series("x", [math.nan, *range(10)], dtype=pl.Float64)
    index = pl.Series("idx", [f"bc{i}" for i in range(len(values))])
    _scale_info, traces = _visualize_numeric(values, index, logger=LOGGER)
    outlier_trace = list(traces)[2]
    assert isinstance(outlier_trace, go.Scattergl)
    assert outlier_trace.y is None or len(outlier_trace.y) == 0


def test_categorical_order_is_sorted() -> None:
    _scale_info, traces = _visualize_categorical(pl.Series("b", [True, None, True, False]))
    (bar_trace,) = traces
    assert isinstance(bar_trace, go.Bar)
    assert bar_trace.x is not None
    assert bar_trace.y is not None
    assert list(bar_trace.x) == [False, True]
    assert list(bar_trace.y) == pytest.approx([100 / 3, 200 / 3])


def test_visualize_numeric_unsigned_and_null() -> None:
    values = pl.Series("x", [None, 0, 1, 1, 2, 2, 3, 10_000], dtype=pl.UInt64)
    index = pl.Series("idx", [f"bc{i}" for i in range(len(values))])
    scale_info, traces = _visualize_numeric(values, index, logger=LOGGER)
    assert scale_info.log_scale is not None
    assert len(list(traces)) == 3

    all_null = pl.Series("x", [None, None], dtype=pl.UInt64)
    scale_info, _traces = _visualize_numeric(all_null, index[:2], logger=LOGGER)
    assert scale_info.unplotable_range


def test_find_qc_files_sorted_and_filtered(tmp_path: Path) -> None:
    for name in ("b.per_cell_qc.tsv.gz", "a.per_cell_qc.tsv.gz", "c.pseudobulk_qc.tsv"):
        (tmp_path / name).touch()
    found = list(_find_qc_files([tmp_path], filter_glob="*per_cell_qc.tsv.gz"))
    assert [path.name for path in found] == [
        "a.per_cell_qc.tsv.gz",
        "b.per_cell_qc.tsv.gz",
    ]


def test_visualize_qc_end_to_end(tmp_path: Path) -> None:
    """Mimic the inputs that WRITE_SUMMARY stages for visualize-qc."""
    table_qc = _write_tsv(
        tmp_path / "metadata.summary-stats.tsv",
        TABLE_COLS,
        [
            ["summary", 2, 1, 1, 0, 1, 0],
            ["IGVFDS0000AAAA", 1, 1, 0, 0, 1, 0],
            ["IGVFDS0000BBBB", 1, 0, 1, 0, 0, 0],
        ],
    )
    accession_dir = tmp_path / "analysis_accession_qc_reports"
    accession_dir.mkdir()
    _write_dict_rows(
        accession_dir / "IGVFDS0000AAAA_per_cell_qc.tsv.gz",
        ACCESSION_COLS_WITH_RNA,
        [_accession_row(f"bc{i}", num_frags=10 * i, rna_read_count=100 * i) for i in range(20)],
    )
    _write_dict_rows(
        accession_dir / "IGVFDS0000BBBB_per_cell_qc.tsv.gz",
        ACCESSION_COLS_ATAC_ONLY,
        [_accession_row(f"bc{i}", num_frags=10 * i, rna_read_count=None) for i in range(20)],
    )
    pseudobulk_dir = tmp_path / "pseudobulk_qc_reports"
    pseudobulk_dir.mkdir()
    accession = "IGVFDS0000AAAA"
    for pseudobulk_id in ("pb1", "pb2"):
        _write_tsv(
            pseudobulk_dir / f"{pseudobulk_id}.per_cell_qc.tsv.gz",
            PSEUDOBULK_COLS,
            [
                # frip is entirely missing, as when CALL_PEAKS did not run for the pseudobulk
                [accession, f"bc{i}", "sub", 100 * i, 50 * i, 1.5, 2.5, i, 0.25, 0.5, 2.0, None]
                for i in range(10)
            ],
        )
    output = tmp_path / "metadata.qc.html"

    visualize_qc(
        output=output,
        table_qc=[table_qc],
        accession_qc=[accession_dir],
        pseudobulk_qc=[pseudobulk_dir],
    )

    html = output.read_text()
    # table, 2 summaries, 2 accessions, 2 pseudobulks
    options = re.findall(r'<option value="fig-(\d+)"', html)
    assert options == [f"{idx}" for idx in range(7)]
    assert html.count('class="chart-container"') == 7
    # plotly.js must be embedded exactly once
    plotly_js = plotly.offline.get_plotlyjs()
    assert html.count(plotly_js) == 1
    figures_html = html.replace(plotly_js, "")
    assert "plotly.js v" not in figures_html
    assert figures_html.count("Plotly.newPlot(") == 6


@pytest.mark.parametrize("title", [None, "explicit"])
def test_add_figure_escapes_title(title: str | None) -> None:
    raw_title = "<b>A & B</b>"
    fig = go.Figure(layout={"title": {"text": raw_title}})
    dropdown_options_html: list[str] = []
    figures_grid_html: list[str] = []
    _add_figure(fig, dropdown_options_html, figures_grid_html, title=title)
    expected_text = "&lt;b&gt;A &amp; B&lt;/b&gt;" if title is None else title
    assert dropdown_options_html == [f'<option value="fig-0" selected>{expected_text}</option>']


def test_visualize_table_css_properties() -> None:
    table_html = _visualize_table(pl.LazyFrame({"a": ["x"], "b": [1]}))
    style = table_html[: table_html.index("</style>")]
    properties = set(re.findall(r"^\s*([\w-]+)\s*:", style, flags=re.MULTILINE))
    assert properties == {"font-weight", "background-color", "border"}
