import dataclasses
import logging
from collections.abc import (
    Iterable,
    Iterator,
)
from enum import Enum
from fnmatch import fnmatch
from pathlib import Path
from typing import (
    ClassVar,
    overload,
)

import numpy as np
import plotly
import plotly.basedatatypes
import plotly.graph_objects as go
import plotly.io as pio
import plotly.subplots as subplots
import polars as pl
from plotly.graph_objs import Figure

# TODO: fix outputs of early pipeline stages to output integer, then change schemas to use integer values
ACCESSION_QC_SCHEMA = pl.Schema(
    {
        "analysis_set_accession": pl.String(),
        "barcode_sample": pl.String(),
        "annotated": pl.Boolean(),
        "found_in_rna": pl.Boolean(),
        "found_in_atac": pl.Boolean(),
        "pseudobulk_id": pl.String(),
        "rna_read_count": pl.Float64(),
        "gene_count": pl.Float64(),
        "pct_mito": pl.Float64(),
        "pct_ribo": pl.Float64(),
        "num_frags": pl.Float64(),
        "pct_duplicated_reads": pl.Float64(),
        "nucleosomal_signal": pl.Float64(),
        "tss_enrichment": pl.Float64(),
    }
)

PSEUDOBULK_QC_SCHEMA = pl.Schema(
    {
        "analysis_set_accession": pl.String(),
        "barcode_sample": pl.String(),
        "subsample": pl.String(),
        "rna_read_count": pl.Float64(),
        "gene_count": pl.Float64(),
        "pct_mito": pl.Float64(),
        "pct_ribo": pl.Float64(),
        "num_frags": pl.Float64(),
        "pct_duplicated_reads": pl.Float64(),
        "nucleosomal_signal": pl.Float64(),
        "tss_enrichment": pl.Float64(),
        "frip": pl.Float64(),
    }
)

TABLE_SCHEMA = pl.Schema(
    {
        "accession": pl.String(),
        "frag": pl.UInt64(),
        "matrix": pl.UInt64(),
        "only-frag": pl.UInt64(),
        "only-matrix": pl.UInt64(),
        "both": pl.UInt64(),
        "neither": pl.UInt64(),
    }
)


@dataclasses.dataclass
class DataScales:
    """Class for storing the scales of a dataset."""

    min_val: int | float
    quartile_1: int | float
    median: int | float
    quartile_3: int | float
    max_val: int | float
    iqr: int | float

    @classmethod
    def from_data(cls, data: pl.Series) -> "DataScales":
        """Produce DataScales from a pandas Series."""
        min_val: float | int = data.min()  # ty: ignore[invalid-assignment]
        quartile_1: float | int
        median: float | int
        quartile_3: float | int
        quartile_1, median, quartile_3 = data.quantile([0.25, 0.5, 0.75])  # ty: ignore[invalid-assignment]
        max_val: float | int = data.max()  # ty: ignore[invalid-assignment]
        return cls(
            min_val=min_val,
            quartile_1=quartile_1,
            median=median,
            quartile_3=quartile_3,
            max_val=max_val,
            iqr=quartile_3 - quartile_1,
        )

    @property
    def high(self) -> int | float:
        """Return the high value of the data range."""
        return min(self.max_val, self.quartile_3 + 1.5 * self.iqr)

    @property
    def low(self) -> int | float:
        """Return the high value of the data range."""
        return max(self.min_val, self.quartile_1 - 1.5 * self.iqr)


def _fmt_text(val: float | int) -> str:
    """Format a value as a string for display in the plot."""
    exp_style = f"{val:.3e}"
    norm_style = f"{int(round(val)):d}" if round(val) == val else f"{val:.3f}"
    return norm_style if len(norm_style) < len(exp_style) else exp_style


@dataclasses.dataclass
class ScaleInfo:
    """Class for storing plotting information of a dataset and managing symlog transformations.

    Attributes:
        min_val: Minimum value of the dataset.
        max_val: Maximum value of the dataset.
        log_scale: Scale where symlog transitions from linear to log.
        num_ticks: Number of ticks for the plot.
        ln_10: Natural logarithm of 10.
        show_x_tick_labels: Whether to show x tick labels.
    """

    min_val: int | float
    max_val: int | float
    log_scale: float | None = None
    num_ticks: int = 6
    ln_10: ClassVar[float] = np.log(10.0)
    show_x_tick_labels: bool = False

    @classmethod
    def from_data_scales(
        cls, data_scales: DataScales, num_ticks: int = 6
    ) -> "ScaleInfo":
        """Produce ScaleInfo from DataScales."""
        if data_scales.max_val - data_scales.quartile_1 > 5 * data_scales.iqr:
            # use a symlog scale, choose the transition from linear to log
            log_scale = max(
                1e-3, min(abs(data_scales.quartile_1), abs(data_scales.quartile_3))
            )
        else:
            log_scale = None
        return cls(
            log_scale=log_scale,
            min_val=data_scales.min_val,
            max_val=data_scales.max_val,
            num_ticks=num_ticks,
        )

    @classmethod
    def _scale(cls, val: int | float, log_scale: float) -> float:
        """Convert unscaled value to symlog. Class method assumes log_scale is not None."""
        abs_val = abs(val)
        return (
            val
            if abs_val < log_scale
            else np.sign(val)
            * (log_scale + np.log1p((abs_val - log_scale) / cls.ln_10) / cls.ln_10)
        )

    @classmethod
    def _inv_scale(cls, val: int | float, log_scale: float) -> float:
        """Convert symlog value to unscaled. Class method assumes log_scale is not None."""
        abs_val = np.abs(val)
        return (
            val
            if abs_val < log_scale
            else np.sign(val)
            * (log_scale + cls.ln_10 * np.expm1(cls.ln_10 * (abs_val - log_scale)))
        )

    def scale(self, val: int | float) -> float:
        """Convert unscaled value to symlog if log_scale is not None."""
        return val if self.log_scale is None else self._scale(val, self.log_scale)

    def inv_scale(self, val: int | float) -> float:
        """Convert symlog value to unscaled if log_scale is not None."""
        return val if self.log_scale is None else self._inv_scale(val, self.log_scale)

    @property
    def type(self) -> str:
        """Return plotly scale type: "-" for linear, "log" for symlog."""
        return "-" if self.log_scale is None else "log"

    @property
    def tickmode(self) -> str:
        """Return tick mode: "auto" for automatic tick placement, "array" for custom tick values."""
        return "auto" if self.tickvals is None else "array"

    @property
    def range(self) -> tuple[float, float]:
        """Return the range of the data: (min_val, max_val)."""
        return self.min_val, self.max_val

    @property
    def transformed_range(self) -> tuple[float, float]:
        """Return the transformed range of the data: (scale(min_val), scale(max_val))."""
        return self.scale(self.min_val), self.scale(self.max_val)

    @property
    def ticktext(self) -> list[str] | None:
        """Return tick text labels for symlog scale, or None for linear scale."""
        if self.log_scale is None:
            return None
        min_scaled_val, max_scaled_val = self.transformed_range
        # get the untransformed tick values that correspond to even divisions in the transformed space,
        # rounded to appropriate precision
        return [
            _fmt_text(self.inv_scale(transformed_tick))
            for transformed_tick in np.linspace(
                min_scaled_val, max_scaled_val, self.num_ticks
            )
        ]

    @property
    def tickvals(self) -> list[float] | None:
        """Return tick values for symlog scale, or None for linear scale."""
        tick_labels = self.ticktext
        return (
            None
            if tick_labels is None
            else [self.scale(float(label)) for label in tick_labels]
        )


def _visualize_table(table_lazy: pl.LazyFrame, title: str) -> str:
    """Create a pure HTML table figure from a DataFrame."""
    table = table_lazy.collect(engine="streaming")

    styled = (
        table.to_pandas()
        .style.hide(axis="index")
        .set_table_attributes(
            'style="border-collapse: collapse; border: 1px solid black;"'
        )
        .set_table_styles(
            [
                # Bold the 1st row (header)
                {
                    "selector": "thead th",
                    "props": "font-weight: bold; font-color: line-color: darkslategray; background-color: Lavender; border: 1px solid black;",
                },
                # Bold the 1st column (index)
                {
                    "selector": "tbody td:first-child",
                    "props": "font-weight: bold; line-color: darkslategray; background-color: Lavender; border: 1px solid black;",
                },
                # General cell styling
                {
                    "selector": "td",
                    "props": "line_color: darkslategray; background-color: lightcyan; border: 1px solid black;",
                },
            ]
        )
    )
    fig = styled.to_html()
    return fig


def _visualize_table_csv(
    path: Path,
    exclude_cols: set[str],
    title: str | None = None,
    schema: pl.Schema | None = None,
) -> go.Figure | str:
    """Read CSV file and create a Plotly table figure from it."""
    table_df = scan_csv(path, exclude_cols=exclude_cols, schema=schema)
    title = path.name if title is None else title
    return _visualize_table(table_df, title=title)


def _get_outliers(
    property_values: pl.Series,
    index: pl.Series,
    low: float,
    high: float,
    scale_info: ScaleInfo,
    max_outliers: int,
) -> tuple[pl.Series, pl.Series]:
    is_outlier = (property_values > high) | (property_values < low)
    outliers = property_values.filter(is_outlier)
    index = index.filter(is_outlier)

    if len(outliers) == 0:
        return outliers, index
    else:
        dtype = outliers.dtype if scale_info.log_scale is None else pl.Float64
        if len(outliers) > max_outliers:
            idx = np.round(np.linspace(0, len(outliers) - 1, max_outliers)).astype(
                np.uint64
            )
            filter_idx = outliers.arg_sort().gather(idx)
            return outliers.gather(filter_idx).map_elements(
                scale_info.scale, return_dtype=dtype
            ), index.gather(filter_idx)
        else:
            return outliers.map_elements(scale_info.scale, return_dtype=dtype), index


def _drop_null(
    property_values: pl.Series, index: pl.Series
) -> tuple[pl.Series, pl.Series]:
    mask = property_values.is_not_null()
    return property_values.filter(mask), index.filter(mask)


def _visualize_numeric(
    property_values: pl.Series, index: pl.Series, max_outliers: int = 5
) -> tuple[ScaleInfo, Iterable[plotly.basedatatypes.BaseTraceType]]:
    """Visualize numeric property values as a box plot with outliers scattered."""
    property_values, index = _drop_null(property_values, index)
    data_scales = DataScales.from_data(property_values)
    scale_info = ScaleInfo.from_data_scales(data_scales)
    high = data_scales.high
    low = data_scales.low

    outliers, outliers_index = _get_outliers(
        property_values,
        index,
        low=low,
        high=high,
        scale_info=scale_info,
        max_outliers=max_outliers,
    )

    unscaled_points = [
        low,
        data_scales.quartile_1,
        data_scales.median,
        data_scales.quartile_3,
        high,
    ]
    scaled_points = [scale_info.scale(_pt) for _pt in unscaled_points]

    box_trace = go.Box(
        x0=0,
        y=[scaled_points[0], scaled_points[-1]],
        name=property_values.name,
        sizemode="quartiles",
        showlegend=False,
        boxpoints=False,
        orientation="v",
        median=[scale_info.scale(data_scales.median)],
        q1=[scale_info.scale(data_scales.quartile_1)],
        q3=[scale_info.scale(data_scales.quartile_3)],
        lowerfence=[scale_info.scale(low)],
        upperfence=[scale_info.scale(high)],
        hoverinfo="skip",
    )
    half_width = 0.25
    scaled_low = scaled_points[0]
    scaled_high = scaled_points[-1]
    # Draw an invisible closed rectangle to display appropriate hover text for the box plot
    # closed rectangle: bottom-left -> bottom-right -> top-right -> top-left -> back
    hover_trace = go.Scatter(
        x=[-half_width, half_width, half_width, -half_width, -half_width],
        y=[scaled_low, scaled_low, scaled_high, scaled_high, scaled_low],
        fill="toself",
        fillcolor="rgba(0,0,0,0)",  # fully transparent fill
        line=dict(color="rgba(0,0,0,0)"),  # invisible border too
        mode="lines",
        name="",  # keep it out of the legend
        showlegend=False,
        text=(
            f"high-whisker: {_fmt_text(unscaled_points[4])}<br>"
            f"quartile 3: {_fmt_text(unscaled_points[3])}<br>"
            f"median: {_fmt_text(unscaled_points[2])}<br>"
            f"quartile 1: {_fmt_text(unscaled_points[1])}<br>"
            f"low-whisker: {_fmt_text(unscaled_points[0])}<br>"
        ),
        hoverinfo="text",  # show only the text on hover
        hoveron="fills",  # trigger hover when over the filled area
    )
    # draw outliers as scatter points
    outlier_trace = go.Scattergl(
        x=[0] * len(outliers),
        y=outliers,
        mode="markers",
        marker=dict(
            size=3,
            color="rgba(90, 90, 90, 0.4)",  # Subtle translucent grey dots
            line=dict(width=0),
        ),
        hoverinfo="text",
        hovertext=outliers_index,
        showlegend=False,
    )

    return scale_info, [box_trace, hover_trace, outlier_trace]


def _visualize_categorical(
    property_values: pl.Series,
) -> tuple[ScaleInfo, Iterable[plotly.basedatatypes.BaseTraceType]]:
    """Visualize categorical property values as a bar chart of percentages."""
    scale_info = ScaleInfo(
        min_val=0.0, max_val=100.0, log_scale=None, show_x_tick_labels=True
    )
    value_counts = property_values.drop_nulls().value_counts()

    bar_trace = go.Bar(
        x=value_counts[property_values.name],
        y=value_counts["count"] * 100 / value_counts["count"].sum(),
        name=property_values.name,
        showlegend=False,
    )

    return scale_info, [bar_trace]


def _visualize_qc(
    qc_df: pl.LazyFrame,
    index: pl.Series,
    title: str,
    logger: logging.Logger,
    max_cols_per_row: int = 6,
    max_outliers: int = 5,
) -> go.Figure:
    """Visualize QC DataFrame as a grid of subplots, one for each column."""

    num_subplots = len(qc_df.collect_schema())
    num_rows = (num_subplots + max_cols_per_row - 1) // max_cols_per_row
    num_cols = min(num_subplots, max_cols_per_row)

    fig = subplots.make_subplots(
        rows=num_rows,
        cols=num_cols,
        shared_yaxes=False,  # Force each column to scale independently
        horizontal_spacing=0.075,  # Tight spacing to maximize horizontal plotting area
        vertical_spacing=0.2,
    )

    for plot_idx, col in enumerate(qc_df.collect_schema().names()):
        logger.debug(f"  plotting {col}")
        property_values = qc_df.select(col).collect(engine="streaming")[col]

        if property_values.dtype.is_numeric():
            scale_info, traces = _visualize_numeric(
                property_values, index, max_outliers=max_outliers
            )
        else:
            scale_info, traces = _visualize_categorical(property_values)

        row = 1 + (plot_idx // num_cols)
        col = 1 + (plot_idx % num_cols)
        for trace in traces:
            fig.add_trace(trace, row=row, col=col)

        fig.update_yaxes(
            nticks=6,  # Limits the density of numeric ticks on the axis
            row=row,
            col=col,
            range=scale_info.transformed_range,
            tickvals=scale_info.tickvals,
            ticktext=scale_info.ticktext,
            tickmode=scale_info.tickmode,
            tickformat=".3g",
        )

        # Hide standard category x-axis labels since subplot titles handle this
        fig.update_xaxes(
            showticklabels=scale_info.show_x_tick_labels,
            row=row,
            col=col,
            title_text=property_values.name,
        )

    # Set layout for whole figure with all subplots
    fig.update_layout(
        title=dict(text=title, font=dict(size=18)),
        template="none",
        height=550,  # Extra vertical height to handle multi-range tracking
        margin=dict(l=100, r=100, t=80, b=40),
        clickmode="event",  # Broadcasts click events out to JavaScript
    )

    return fig


def _visualize_qc_csv(
    path: Path,
    exclude_cols: set[str],
    index_col: str,
    logger: logging.Logger,
    title: str | None = None,
    schema: pl.Schema | None = None,
    max_outliers: int = 5,
    max_cols_per_row: int = 6,
) -> go.Figure:
    """Read CSV file and visualize QC DataFrame as a grid of subplots."""
    qc_df, index = scan_csv(
        path, exclude_cols=exclude_cols, index_col=index_col, schema=schema
    )
    title = path.name.split(".", 1)[0] if title is None else title
    return _visualize_qc(
        qc_df=qc_df,
        index=index,
        title=title,
        logger=logger,
        max_outliers=max_outliers,
        max_cols_per_row=max_cols_per_row,
    )


def _visualize_summary_qc(
    paths: Iterable[Path],
    title: str,
    exclude_cols: set[str],
    index_col: str,
    logger: logging.Logger,
    schema: pl.Schema | None = None,
    max_outliers: int = 5,
    max_cols_per_row: int = 6,
) -> go.Figure:
    """Visualize summary QC by concatenating dataframes from multiple CSV files."""
    lazy_dfs, indices = zip(
        *(
            scan_csv(
                path, exclude_cols=exclude_cols, index_col=index_col, schema=schema
            )
            for path in paths
        )
    )
    qc_df = pl.concat(lazy_dfs, how="diagonal_relaxed")
    index = pl.concat(indices, how="vertical")
    return _visualize_qc(
        qc_df=qc_df,
        index=index,
        title=title,
        logger=logger,
        max_outliers=max_outliers,
        max_cols_per_row=max_cols_per_row,
    )


@overload
def scan_csv(
    path: Path,
    index_col: str,
    sep: str | None = None,
    exclude_cols: set[str] | None = None,
    schema: pl.Schema | None = None,
) -> tuple[pl.LazyFrame, pl.Series]: ...


@overload
def scan_csv(
    path: Path,
    index_col: None = None,
    sep: str | None = None,
    exclude_cols: set[str] | None = None,
    schema: pl.Schema | None = None,
) -> pl.LazyFrame: ...


def scan_csv(
    path: Path,
    index_col: str | None = None,
    sep: str | None = None,
    exclude_cols: set[str] | None = None,
    schema: pl.Schema | None = None,
) -> pl.LazyFrame | tuple[pl.LazyFrame, pl.Series]:
    """Read CSV file with optional separator inference and keyword arguments."""
    if sep is None:
        match path.suffixes:
            case [*_parts, ".tsv"] | [*_parts, ".tsv", ".gz"]:
                sep = "\t"
            case [*_parts, ".csv"] | [*_parts, ".csv", ".gz"]:
                sep = ","
            case _:
                raise ValueError(
                    f"Could not infer separator for file {path} with suffix {path.suffix}. Please "
                    "specify separator with 'sep' argument."
                )
    df_lazy = pl.scan_csv(path, separator=sep, infer_schema=False, schema=schema)
    schema = df_lazy.collect_schema()
    if index_col is None:
        return (
            df_lazy
            if exclude_cols is None
            else df_lazy.drop(exclude_cols.intersection(schema))
        )
    else:
        if index_col not in schema:
            raise ValueError(f"Invalid index column: {index_col}")
        index = df_lazy.select(index_col).collect(engine="streaming")[index_col]
        if exclude_cols is None:
            exclude_cols = {index_col}
        else:
            exclude_cols.add(index_col)
        df_lazy = df_lazy.drop(exclude_cols.intersection(schema))
        return df_lazy, index


def _is_csv_or_tsv(path: Path) -> bool:
    """Check if the file is a CSV or TSV file."""
    return path.is_file() and path.suffix in [".csv", ".tsv", ".csv.gz", ".tsv.gz"]


def _find_qc_files(input_paths: list[Path], filter_glob: str = "*") -> Iterator[Path]:
    """Yield paths of CSV/TSV files matching the filter glob from the input paths."""
    for input_path in input_paths:
        if input_path.is_file():
            if fnmatch(f"{input_path}", filter_glob):
                yield input_path
        else:
            for filtered_path in input_path.glob(filter_glob):
                if _is_csv_or_tsv(filtered_path):
                    yield filtered_path


def _get_master_template(
    dropdown_options_html: list[str], figures_grid_html: list[str]
) -> str:
    """Return the master HTML template for the QC dashboard."""
    # Double the outer brackets to properly escape them inside the Python f-string
    return f"""<!DOCTYPE html>
    <html>
    <head>
        <title>QC Dashboard</title>
        <script src="https://plot.ly"></script>
        <style>
            body {{
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                margin: 20px;
                background-color: #f8f9fa;
            }}
            .header-controls {{
                background: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
                margin-bottom: 20px;
            }}
            select {{
                padding: 10px 15px;
                font-size: 16px;
                border-radius: 5px;
                border: 1px solid #ccc;
                width: 350px;
            }}
            .dashboard-grid {{
                display: grid;
                grid-template-columns: 1fr;
                gap: 20px;
            }}
            .chart-container {{
                background: white;
                padding: 15px;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
                position: relative;
            }}
            #toast-notification {{
                position: fixed;
                top: 20px;
                right: 20px;
                background: #28a745;
                color: white;
                padding: 12px 22px;
                border-radius: 6px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.15);
                font-size: 14px;
                font-weight: 500;
                opacity: 0;
                transition: opacity 0.2s ease;
                pointer-events: none;
                z-index: 10000;
            }}
            .js-plotly-plot .plotly .cursor-pointer {{
                cursor: copy !important;
            }}
        </style>
    </head>
    <body>

        <div id="toast-notification">Copied to Clipboard!</div>

        <div class="header-controls">
            <label for="figure-selector" style="font-weight: bold; margin-right: 10px;">Select Visual Analysis:</label>
            <select id="figure-selector" onchange="switchFigure()">
                {"".join(dropdown_options_html)}
            </select>
        </div>

        <div class="dashboard-grid">
            {"".join(figures_grid_html)}
        </div>

        <script>
        function switchFigure() {{
            const selector = document.getElementById('figure-selector');
            const targetId = selector.value;
            const allCharts = document.querySelectorAll('.chart-container');

            allCharts.forEach(chart => {{
                if (chart.id === targetId) {{
                    chart.style.display = 'block';
                    window.dispatchEvent(new Event('resize'));
                }} else {{
                    chart.style.display = 'none';
                }}
            }});
        }}

        document.addEventListener('DOMContentLoaded', function () {{
            function copyTextToClipboard(text) {{
                if (navigator.clipboard && navigator.clipboard.writeText) {{
                    navigator.clipboard.writeText(text).then(showToast, () => fallbackCopy(text));
                }} else {{
                    fallbackCopy(text);
                }}
            }}

            function fallbackCopy(text) {{
                const ta = document.createElement('textarea');
                ta.value = text;
                ta.style.position = 'fixed';
                document.body.appendChild(ta);
                ta.select();
                try {{
                    document.execCommand('copy');
                    showToast();
                }} catch (err) {{
                    console.error('Fallback copy failed', err);
                }}
                document.body.removeChild(ta);
            }}

            function showToast() {{
                const toast = document.getElementById('toast-notification');
                toast.style.opacity = '1';
                setTimeout(() => {{ toast.style.opacity = '0'; }}, 1500);
            }}

            const graphDivs = document.querySelectorAll('.plotly-graph-div');
            graphDivs.forEach(function (graphDiv) {{
                graphDiv.on('plotly_click', function (data) {{
                    if (!data || !data.points || data.points.length === 0) return;
                    const point = data.points[0];
                    if (point.text !== undefined && point.text !== '') {{
                        copyTextToClipboard(point.text);
                    }} else if (point.hovertext !== undefined && point.hovertext !== '') {{
                        copyTextToClipboard(point.hovertext);
                    }}
                }});
            }});
        }});
        </script>
    </body>
    </html>"""


# Ensure CDN is used to extract clean, minimal script snippets without bloated engine duplication
pio.renderers.default = "browser"


def _add_figure(
    fig: Figure | str,
    dropdown_options_html: list[str],
    figures_grid_html: list[str],
    title: str | None = None,
) -> None:
    """Add a figure to the dropdown options and figures grid HTML lists."""
    # 2. Extract ONLY the necessary raw <div> components of this individual chart
    # setting full_html=False drops the repetitive <html> headers, keeping file footprint tiny
    fig_div = (
        fig
        if isinstance(fig, str)
        else pio.to_html(fig, full_html=False, include_plotlyjs="none")
    )
    title_text = (
        title
        if title is not None
        else fig.layout.title.text
        if isinstance(fig, Figure)
        else "Summary stats"
    )

    # 3. Add to dropdown options (Value points to our specific element ID)
    # The first item will be marked selected
    index = len(dropdown_options_html)
    selected_attr = "selected" if index == 0 else ""
    dropdown_options_html.append(
        f'<option value="fig-{index}" {selected_attr}>{title_text}</option>'
    )

    # 4. Wrap the graph in an isolated, queryable HTML component container
    # Hide all figures by default except for the very first one
    display_style = "block" if index == 0 else "none"
    figures_grid_html.append(
        f'<div id="fig-{index}" class="chart-container" style="display: {display_style};">{fig_div}</div>'
    )


def _find_input_qc_files(
    table_qc: list[Path],
    accession_qc: list[Path],
    pseudobulk_qc: list[Path],
    logger: logging.Logger,
) -> tuple[list[Path], list[Path], list[Path]]:
    """Find the paths to all the CSV/TSVs that will be visualized."""
    table_qc_files: list[Path] = list(
        _find_qc_files(input_paths=table_qc, filter_glob="*.summary-stats.tsv")
    )
    logger.info(f"Found {len(table_qc_files)} table QC files:")
    for file in table_qc_files:
        logger.info(f"  {file}")

    accession_qc_files: list[Path] = list(
        _find_qc_files(input_paths=accession_qc, filter_glob="*per_cell_qc.tsv")
    )
    logger.info(f"Found {len(accession_qc_files)} accession QC files:")
    for file in accession_qc_files:
        logger.info(f"  {file}")

    pseudobulk_qc_files: list[Path] = list(
        _find_qc_files(input_paths=pseudobulk_qc, filter_glob="*per_cell_qc.tsv")
    )
    logger.info(f"Found {len(pseudobulk_qc_files)} pseudobulk QC files:")
    for file in pseudobulk_qc_files:
        logger.info(f"  {file}")
    return table_qc_files, accession_qc_files, pseudobulk_qc_files


class LogLevel(Enum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


def visualize_qc(
    *,
    output: Path,
    table_qc: list[Path] = [],
    accession_qc: list[Path] = [],
    pseudobulk_qc: list[Path] = [],
    exclude_col: list[str] = ["analysis_set_accession", "pseudobulk_id", "subsample"],
    index_col: str = "barcode_sample",
    max_outliers: int = 5,
    max_cols_per_row: int = 6,
    plot_accessions_summary: bool | None = None,
    plot_pseudobulks_summary: bool | None = None,
    plot_individual_accessions: bool | None = None,
    plot_individual_pseudobulks: bool | None = None,
    log_level: LogLevel = LogLevel.INFO,
) -> None:
    """
    Visualize QC metrics for inputs summary table, accession, and pseudobulk QC files.

    Writes an interactive HTML report to the output path.

    Args:
        output: Path to the output HTML file
        table_qc: One or more paths to table QC files
        accession_qc: One or more paths to accession QC files
        pseudobulk_qc: One or more paths to pseudobulk QC files
        exclude_col: Comma-separatedlist of column names to exclude from visualizing in the QC plots
        index_col: Name of the index column (will show up as hover text for outliers in the QC plots)
        max_outliers: Maximum number of outliers to display in the QC plots. More results in more memory usage.
        max_cols_per_row: Maximum number of columns per row in the QC plots.
        plot_accessions_summary: Whether to plot a summary of accession QC metrics.
            If None, will be set to True if there are multiple accession QC files, False otherwise.
        plot_pseudobulks_summary: Whether to plot a summary of pseudobulk QC metrics.
            If None, will be set to True if there are multiple pseudobulk QC files, False otherwise.
        plot_individual_accessions: Whether to plot individual accession QC metrics.
            If None, will be set to True if there are less than 10 accession QC files, False otherwise.
        plot_individual_pseudobulks: Whether to plot individual pseudobulk QC metrics.
            If None, will be set to True if there are less than 10 pseudobulk QC files, False otherwise.
        log_level: Logging level to report.
    """
    logging.basicConfig(
        level=log_level.value,
        format="%(asctime)s %(name)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger("visualize-qc")
    table_qc_files, accession_qc_files, pseudobulk_qc_files = _find_input_qc_files(
        table_qc=table_qc,
        accession_qc=accession_qc,
        pseudobulk_qc=pseudobulk_qc,
        logger=logger,
    )

    if plot_accessions_summary is None:
        plot_accessions_summary = len(accession_qc_files) > 1
    if plot_pseudobulks_summary is None:
        plot_pseudobulks_summary = len(pseudobulk_qc_files) > 1
    if plot_individual_accessions is None:
        plot_individual_accessions = len(accession_qc_files) < 10
    if plot_individual_pseudobulks is None:
        plot_individual_pseudobulks = len(pseudobulk_qc_files) < 10

    exclude_cols = {_split_col for col in exclude_col for _split_col in col.split(",")}

    # Storage lists for raw HTML parts
    dropdown_options_html: list[str] = []
    figures_grid_html: list[str] = []

    for table_qc_file in table_qc_files:
        logger.info(f"Visualizing table QC file: {table_qc_file}")
        _fig = _visualize_table_csv(
            table_qc_file,
            exclude_cols=exclude_cols,
            schema=TABLE_SCHEMA,
        )
        _add_figure(
            _fig,
            dropdown_options_html,
            figures_grid_html,
            title=table_qc_file.name.split(".", 1)[0],
        )

    if len(accession_qc_files) > 0 and plot_accessions_summary:
        logger.info("Visualizing summary accession QC")
        _fig = _visualize_summary_qc(
            accession_qc_files,
            exclude_cols=exclude_cols,
            index_col=index_col,
            title="All accession QCs",
            logger=logger,
            schema=ACCESSION_QC_SCHEMA,
            max_outliers=max_outliers,
            max_cols_per_row=max_cols_per_row,
        )
        _add_figure(_fig, dropdown_options_html, figures_grid_html)
    if len(pseudobulk_qc_files) > 0 and plot_pseudobulks_summary:
        logger.info("Visualizing summary pseudobulk QC")
        _fig = _visualize_summary_qc(
            pseudobulk_qc_files,
            exclude_cols=exclude_cols,
            index_col=index_col,
            title="All pseudobulk QCs",
            logger=logger,
            schema=PSEUDOBULK_QC_SCHEMA,
            max_outliers=max_outliers,
            max_cols_per_row=max_cols_per_row,
        )
        _add_figure(_fig, dropdown_options_html, figures_grid_html)

    if plot_individual_accessions:
        for idx, input_qc_file in enumerate(accession_qc_files):
            logger.info(
                f"Visualizing accession QC file: {input_qc_file} ({idx + 1}/{len(accession_qc_files)})"
            )
            _fig = _visualize_qc_csv(
                input_qc_file,
                exclude_cols=exclude_cols,
                index_col=index_col,
                logger=logger,
                schema=ACCESSION_QC_SCHEMA,
                max_outliers=max_outliers,
                max_cols_per_row=max_cols_per_row,
            )
            _add_figure(_fig, dropdown_options_html, figures_grid_html)

    if plot_individual_pseudobulks:
        for idx, input_qc_file in enumerate(pseudobulk_qc_files):
            logger.info(
                f"Visualizing pseudobulk QC file: {input_qc_file} ({idx + 1}/{len(pseudobulk_qc_files)})"
            )
            _fig = _visualize_qc_csv(
                input_qc_file,
                exclude_cols=exclude_cols,
                index_col=index_col,
                title=input_qc_file.parent.name,
                logger=logger,
                schema=PSEUDOBULK_QC_SCHEMA,
                max_outliers=max_outliers,
                max_cols_per_row=max_cols_per_row,
            )
            _add_figure(_fig, dropdown_options_html, figures_grid_html)

    # Create the Master HTML Dashboard Framework
    master_html_template = _get_master_template(
        dropdown_options_html=dropdown_options_html, figures_grid_html=figures_grid_html
    )

    # Write the unified multi-chart wrapper to HTML
    with output.open("w", encoding="utf-8") as f:
        f.write(master_html_template)
    logger.info(f"Output written to: {output}")
