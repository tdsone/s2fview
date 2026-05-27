"""Plotly-based alternative to the matplotlib viewer in :mod:`s2fview`.

Why this exists
---------------
The matplotlib + ipympl path renders every frame Python-side and pushes a
PNG diff over the WebSocket. Hover, zoom, and pan therefore all incur a
round-trip per interaction — fine for casual use but never going to feel
truly instant on a remote kernel.

Plotly renders client-side in WebGL/Canvas. Hover crosshairs (spike lines),
zoom, pan, and reset are all native browser operations — no Python round-
trip. The trade-off: less fine-grained styling control than matplotlib, a
heavier JS bundle on first display, and the viewer doesn't share matplotlib
artists with the static plots.

API
---
:func:`plotly_coverage_track` mirrors :func:`s2fview.coverage_track` as
closely as possible (same `Gene` dataclass, same kwargs where meaningful)
so you can A/B by swapping the import.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

from . import SEQUENCE_COLORS, Gene, _assign_lanes

if TYPE_CHECKING:
    import plotly.graph_objects as go  # noqa: F401


_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
_SEQ_COLORS_ORDERED = [
    SEQUENCE_COLORS["A"],
    SEQUENCE_COLORS["C"],
    SEQUENCE_COLORS["G"],
    SEQUENCE_COLORS["T"],
    SEQUENCE_COLORS["N"],
]


def plotly_coverage_track(
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
    """Plotly equivalent of :func:`s2fview.coverage_track`.

    Renders entirely client-side; hover, pan, zoom, and reset are native
    Plotly interactions with no Python round-trip. Returns a
    ``plotly.graph_objects.Figure``.
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
        row_heights.append(0.12)
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
        fig.update_yaxes(
            showticklabels=False, showgrid=False, zeroline=False, row=seq_row, col=1
        )

    if genes:
        gene_row = n_rows
        _add_gene_shapes(
            fig,
            genes,
            lanes,
            n_lanes,
            forward_color,
            reverse_color,
            gene_row,
            x_range=(float(positions[0]), float(positions[-1])),
            show_chevrons=show_chevrons,
        )
        fig.update_yaxes(
            range=[n_lanes, 0],  # lane 0 at top, like matplotlib version
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            row=gene_row,
            col=1,
        )

    # Crosshair: spike lines on the shared x-axis, drawn natively by Plotly
    # in the browser. No Python round-trip, no lag.
    fig.update_xaxes(
        showspikes=True,
        spikemode="across",
        spikethickness=1,
        spikedash="dot",
        spikecolor=crosshair_color,
        showline=False,
    )

    # Tick labels only on the coverage axis (top); other rows hide them.
    fig.update_xaxes(showticklabels=True, row=1, col=1)
    for r in range(2, n_rows + 1):
        fig.update_xaxes(showticklabels=False, row=r, col=1)

    if height is None:
        height = 300 + (80 if sequence is not None else 0) + (40 * max(1, n_lanes) if genes else 0)

    fig.update_layout(
        title=title,
        width=width,
        height=height,
        hovermode="x",
        showlegend=False,
        margin={"l": 60, "r": 20, "t": 50 if title else 20, "b": 40},
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": "Inter, -apple-system, Helvetica, Arial, sans-serif"},
        modebar={"orientation": "v"},
    )

    fig.update_yaxes(title_text="Coverage", row=1, col=1, zeroline=True, zerolinecolor="#cccccc")

    return fig


def _add_gene_shapes(
    fig,
    genes: Sequence[Gene],
    lanes: Sequence[int],
    n_lanes: int,
    forward_color: str,
    reverse_color: str,
    row: int,
    x_range: tuple[float, float],
    show_chevrons: bool,
) -> None:
    """Add gene boxes (with directional arrow tips) + chevrons to a subplot row.

    Each gene is one ``Scatter`` trace with ``fill='toself'`` so the gene
    shapes hover-test as a unit and we get a nice per-gene tooltip.
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

        # Thin intron line (full span) as a scatter line trace.
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

        # Exon polygons -- terminal exon (in transcription direction) becomes
        # a pentagon with a directional arrow tip.
        exons = sorted(gene.exons) if gene.exons else [(gene.start, gene.end)]
        terminal_idx = len(exons) - 1 if gene.strand == "+" else 0
        xs: list[float] = []
        ys: list[float] = []
        for i, (ex_start, ex_end) in enumerate(exons):
            if i > 0:
                xs.append(None)  # disconnect from previous exon polygon
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

        # Gene name label
        fig.add_annotation(
            x=(gene.start + gene.end) / 2,
            y=lane + 0.12,
            text=gene.name,
            showarrow=False,
            font={"color": gene_color, "size": 10},
            xref=f"x{row}" if row > 1 else "x",
            yref=f"y{row}" if row > 1 else "y",
        )

        # Strand chevrons along intron line, drawn under the exon boxes by
        # rendering them as a separate fill='toself' trace ordered first.
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
    """Return polygon vertices for an exon (rectangle or pentagon)."""
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
    # left
    body_start = start + arrow_width
    return [
        (body_start, y0),
        (end, y0),
        (end, y1),
        (body_start, y1),
        (start, y_mid),
        (body_start, y0),
    ]
