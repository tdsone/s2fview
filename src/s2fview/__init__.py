"""s2fview: a Plotly-based genome track viewer for Jupyter.

The viewer renders coverage + DNA sequence + gene annotations as a stack of
client-side Plotly subplots. Hover (spike-line crosshair), pan, zoom, and
reset are all native browser interactions — there is no Python round-trip
per mouse move, so it stays smooth on remote kernels.

The main entrypoint is :func:`coverage_track`. :class:`Gene` is a tiny
dataclass for gene annotations.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import plotly.graph_objects as go  # noqa: F401


# IGV-ish palette for DNA bases. Saturated enough to be unambiguous, soft
# enough that white letters on top read cleanly.
SEQUENCE_COLORS: dict[str, str] = {
    "A": "#3aa55a",  # green
    "C": "#3585c7",  # blue
    "G": "#e8a43c",  # orange
    "T": "#d65a4a",  # red
    "N": "#a8a8a8",  # gray fallback
}

_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
_SEQ_COLORS_ORDERED = [
    SEQUENCE_COLORS["A"],
    SEQUENCE_COLORS["C"],
    SEQUENCE_COLORS["G"],
    SEQUENCE_COLORS["T"],
    SEQUENCE_COLORS["N"],
]

_FONT_STACK = "Inter, -apple-system, Helvetica, Arial, sans-serif"


@dataclass(frozen=True)
class Gene:
    """A minimal gene annotation.

    Parameters
    ----------
    name, start, end, strand:
        Self-explanatory.
    exons:
        Optional list of ``(exon_start, exon_end)`` pairs. When provided, the
        gene is drawn as thick exon boxes joined by a thin intron line spanning
        ``start..end``. When ``None``, the gene is rendered as a single solid
        box from ``start`` to ``end``.
    """

    name: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"
    exons: tuple[tuple[int, int], ...] | None = field(default=None)


def coverage_track(
    values: Sequence[float],
    positions: Sequence[int] | None = None,
    *,
    genes: Sequence[Gene] | None = None,
    sequence: str | None = None,
    color: str = "#3b7dd8",
    forward_color: str = "#2563eb",
    reverse_color: str = "#dc2626",
    crosshair_color: str = "#555555",
    title: str | None = None,
    width: int = 950,
    height: int | None = None,
    show_chevrons: bool = True,
):
    """Build a Plotly viewer with coverage + optional sequence + gene tracks.

    All interactions run client-side — hover spike-lines, pan, zoom, and the
    standard Plotly modebar. The figure is a stack of x-shared subplots:

    * Coverage (top) — filled step area, with tick labels.
    * DNA sequence (middle, when ``sequence`` is provided) — colored cells
      per base, with bold letters that appear on zoom-in.
    * Gene track (bottom, when ``genes`` is provided) — exon boxes with a
      directional arrow tip on the terminal exon, thin intron lines with
      strand chevrons, and gene-name labels color-coded by strand.

    Returns either a ``plotly.graph_objects.Figure`` or, when ``sequence``
    is provided, a ``plotly.graph_objects.FigureWidget`` (needed for the
    on-zoom letter-visibility callback). Both render directly in Jupyter
    via their ``_repr_mimebundle_``.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if positions is None:
        positions = list(range(len(values)))
    if len(positions) != len(values):
        raise ValueError(
            f"positions and values must have the same length "
            f"(got {len(positions)} and {len(values)})"
        )
    if sequence is not None and len(sequence) != len(values):
        raise ValueError(
            f"sequence length must match values length "
            f"(got {len(sequence)} and {len(values)})"
        )

    rows: list[str] = ["coverage"]
    row_heights: list[float] = [1.0]
    if sequence is not None:
        rows.append("sequence")
        row_heights.append(0.25)
    lanes: list[int] = []
    n_lanes = 0
    if genes:
        lanes = _assign_lanes(genes)
        n_lanes = max(lanes) + 1
        rows.append("genes")
        row_heights.append(0.10 * n_lanes + 0.05)

    n_rows = len(rows)
    fig = make_subplots(
        rows=n_rows,
        cols=1,
        shared_xaxes=True,
        row_heights=row_heights,
        vertical_spacing=0.03,
    )

    fig.add_trace(
        go.Scatter(
            x=list(positions),
            y=list(values),
            mode="lines",
            line={"color": color, "width": 1, "shape": "hv"},
            fill="tozeroy",
            fillcolor=color,
            opacity=0.85,
            name="coverage",
            hovertemplate="x = %{x}<br>coverage = %{y}<extra></extra>",
        ),
        row=1,
        col=1,
    )

    if sequence is not None:
        seq_row = 2
        indices = [[_BASE_TO_IDX.get(b.upper(), 4) for b in sequence]]
        # Discrete colorscale so each integer maps to one color.
        n = len(_SEQ_COLORS_ORDERED)
        colorscale = []
        for i, c in enumerate(_SEQ_COLORS_ORDERED):
            colorscale.append([i / n, c])
            colorscale.append([(i + 1) / n, c])
        fig.add_trace(
            go.Heatmap(
                z=indices,
                x=list(positions),
                zmin=-0.5,
                zmax=n - 0.5,
                colorscale=colorscale,
                showscale=False,
                hovertemplate="x = %{x}<br>base = %{customdata}<extra></extra>",
                customdata=[list(sequence)],
            ),
            row=seq_row,
            col=1,
        )
        # Letter glyphs as a Scatter text trace on top of the heatmap. Hidden
        # initially; the FigureWidget callback below flips visibility based on
        # visible base density.
        seq_upper = sequence.upper()
        fig.add_trace(
            go.Scatter(
                x=list(positions),
                y=[0.0] * len(positions),
                mode="text",
                text=list(seq_upper),
                textposition="middle center",
                textfont={
                    "family": _FONT_STACK,
                    "size": 13,
                    "color": "white",
                    "weight": "bold",
                },
                hoverinfo="skip",
                showlegend=False,
                name="sequence_letters",
            ),
            row=seq_row,
            col=1,
        )
        fig.update_yaxes(
            range=[-0.5, 0.5],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=seq_row,
            col=1,
        )

    if genes:
        gene_row = n_rows
        _add_gene_shapes(
            fig,
            genes,
            lanes,
            forward_color,
            reverse_color,
            gene_row,
            x_range=(float(positions[0]), float(positions[-1])),
            show_chevrons=show_chevrons,
        )
        fig.update_yaxes(
            range=[n_lanes, 0],  # lane 0 at top
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=gene_row,
            col=1,
        )

    # Native spike-line crosshair on the shared x-axis — drawn by Plotly
    # client-side, no Python round-trip.
    fig.update_xaxes(
        showspikes=True,
        spikemode="across",
        spikethickness=1,
        spikedash="dot",
        spikecolor=crosshair_color,
        showline=False,
    )

    # Tick labels on the coverage axis only.
    fig.update_xaxes(showticklabels=True, row=1, col=1)
    for r in range(2, n_rows + 1):
        fig.update_xaxes(showticklabels=False, row=r, col=1)

    if height is None:
        height = (
            300
            + (80 if sequence is not None else 0)
            + (40 * max(1, n_lanes) if genes else 0)
        )

    fig.update_layout(
        title=title,
        width=width,
        height=height,
        hovermode="x",
        showlegend=False,
        margin={"l": 60, "r": 20, "t": 50 if title else 20, "b": 40},
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": _FONT_STACK},
        modebar={"orientation": "v"},
    )

    fig.update_yaxes(title_text="Coverage", row=1, col=1, zeroline=True, zerolinecolor="#cccccc")

    if sequence is None:
        return fig

    # Lazy letters: wrap as a FigureWidget and bind a relayout callback that
    # toggles the letter trace's visibility based on visible base density.
    fw = go.FigureWidget(fig)
    letter_trace_idx: int | None = None
    for i, tr in enumerate(fw.data):
        if getattr(tr, "name", None) == "sequence_letters":
            letter_trace_idx = i
            break

    if letter_trace_idx is None:
        return fw

    plot_width_px = max(1, width - 80)
    visible_bases_threshold = max(8, plot_width_px // 10)  # ~10 px per letter

    def _toggle_letters(layout, x_range=None):
        if x_range is None:
            xr = getattr(fw.layout.xaxis, "range", None)
            if xr is None:
                xr = (float(positions[0]) - 0.5, float(positions[-1]) + 0.5)
            x_range = (float(xr[0]), float(xr[1]))
        visible = (x_range[1] - x_range[0]) <= visible_bases_threshold
        with fw.batch_update():
            fw.data[letter_trace_idx].visible = bool(visible)

    fw.data[letter_trace_idx].visible = False
    fw.layout.on_change(_toggle_letters, "xaxis.range")

    return fw


def _add_gene_shapes(
    fig,
    genes: Sequence[Gene],
    lanes: Sequence[int],
    forward_color: str,
    reverse_color: str,
    row: int,
    x_range: tuple[float, float],
    show_chevrons: bool,
) -> None:
    """Draw gene boxes (with arrow tips) + chevrons + labels in ``row``.

    Each gene gets one ``Scatter`` polygon trace for its exons (so the box
    hover-tests as a unit with a per-gene tooltip), a thin line trace for
    introns, and an annotation for the name.
    """
    import plotly.graph_objects as go

    span = x_range[1] - x_range[0]
    arrow_w = max(1.0, span * 0.012)
    chev_half_w = max(0.5, span * 0.005)
    chev_half_h = 0.12
    box_h = 0.55

    for gene, lane in zip(genes, lanes):
        gene_color = forward_color if gene.strand == "+" else reverse_color
        center_y = lane + 0.5
        y0 = center_y - box_h / 2
        y1 = center_y + box_h / 2

        fig.add_trace(
            go.Scatter(
                x=[gene.start, gene.end],
                y=[center_y, center_y],
                mode="lines",
                line={"color": gene_color, "width": 1.5},
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row,
            col=1,
        )

        exons = sorted(gene.exons) if gene.exons else [(gene.start, gene.end)]
        terminal_idx = len(exons) - 1 if gene.strand == "+" else 0
        xs: list[float] = []
        ys: list[float] = []
        for i, (ex_start, ex_end) in enumerate(exons):
            if i > 0:
                xs.append(None)
                ys.append(None)
            verts = _exon_vertices(
                ex_start,
                ex_end,
                y0,
                y1,
                center_y,
                arrow_side=("right" if gene.strand == "+" else "left")
                if i == terminal_idx
                else None,
                arrow_width=min(arrow_w, (ex_end - ex_start) * 0.5),
            )
            for vx, vy in verts:
                xs.append(vx)
                ys.append(vy)

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                fill="toself",
                fillcolor=gene_color,
                mode="lines",
                line={"color": gene_color, "width": 0},
                hovertemplate=(
                    f"<b>{gene.name}</b><br>"
                    f"{gene.start:,}–{gene.end:,} ({gene.strand})<extra></extra>"
                ),
                showlegend=False,
            ),
            row=row,
            col=1,
        )

        fig.add_annotation(
            x=(gene.start + gene.end) / 2,
            y=lane + 0.12,
            text=gene.name,
            showarrow=False,
            font={"color": gene_color, "size": 10},
            xref=f"x{row}" if row > 1 else "x",
            yref=f"y{row}" if row > 1 else "y",
        )

        if show_chevrons:
            chev_xs: list[float] = []
            chev_ys: list[float] = []
            span_g = gene.end - gene.start
            n_chev = max(1, int(span_g / 35))
            for k in range(n_chev):
                cx = gene.start + (k + 0.5) * span_g / n_chev
                if gene.strand == "+":
                    tri = [
                        (cx - chev_half_w, center_y - chev_half_h),
                        (cx + chev_half_w, center_y),
                        (cx - chev_half_w, center_y + chev_half_h),
                    ]
                else:
                    tri = [
                        (cx + chev_half_w, center_y - chev_half_h),
                        (cx - chev_half_w, center_y),
                        (cx + chev_half_w, center_y + chev_half_h),
                    ]
                if chev_xs:
                    chev_xs.append(None)
                    chev_ys.append(None)
                for vx, vy in tri:
                    chev_xs.append(vx)
                    chev_ys.append(vy)
            if chev_xs:
                fig.add_trace(
                    go.Scatter(
                        x=chev_xs,
                        y=chev_ys,
                        fill="toself",
                        fillcolor=gene_color,
                        mode="lines",
                        line={"color": gene_color, "width": 0},
                        hoverinfo="skip",
                        showlegend=False,
                    ),
                    row=row,
                    col=1,
                )


def _exon_vertices(
    start: float,
    end: float,
    y0: float,
    y1: float,
    y_mid: float,
    *,
    arrow_side: str | None,
    arrow_width: float,
) -> list[tuple[float, float]]:
    """Polygon vertices for an exon — rectangle, or pentagon with an arrow tip."""
    if not arrow_side or arrow_width <= 0:
        return [(start, y0), (end, y0), (end, y1), (start, y1), (start, y0)]
    arrow_width = min(arrow_width, end - start)
    if arrow_side == "right":
        body_end = end - arrow_width
        return [
            (start, y0),
            (body_end, y0),
            (end, y_mid),
            (body_end, y1),
            (start, y1),
            (start, y0),
        ]
    body_start = start + arrow_width
    return [
        (body_start, y0),
        (end, y0),
        (end, y1),
        (body_start, y1),
        (start, y_mid),
        (body_start, y0),
    ]


def _assign_lanes(genes: Sequence[Gene]) -> list[int]:
    """Greedy interval-packing: assign each gene to the lowest free lane."""
    order = sorted(range(len(genes)), key=lambda i: (genes[i].start, genes[i].end))
    lane_ends: list[int] = []
    result = [0] * len(genes)
    for i in order:
        g = genes[i]
        placed = False
        for lane, end in enumerate(lane_ends):
            if g.start >= end:
                lane_ends[lane] = g.end
                result[i] = lane
                placed = True
                break
        if not placed:
            result[i] = len(lane_ends)
            lane_ends.append(g.end)
    return result
