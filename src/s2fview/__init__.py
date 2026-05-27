from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from matplotlib.widgets import SpanSelector

# Prefer a nicer sans-serif font chain. Matplotlib walks the list and uses the
# first installed one, silently falling back to the platform default otherwise.
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = [
    "Avenir Next",
    "Inter",
    "SF Pro Text",
    "Helvetica Neue",
    "Helvetica",
    *mpl.rcParams["font.sans-serif"],
]


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
    color: str = "#3b7dd8",
    forward_color: str = "#2563eb",
    reverse_color: str = "#dc2626",
    label: str | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str = "Coverage",
    figsize: tuple[float, float] = (10, 2.2),
    dpi: int = 200,
    guide_color: str = "#9a9a9a",
) -> tuple[Figure, Axes]:
    """Plot a coverage track, optionally with a gene track below the x-axis.

    When ``genes`` is provided, the figure has two stacked, x-shared axes:
    the coverage on top (keeping its tick labels), and a dedicated gene-
    annotation strip beneath. Overlapping genes are automatically stacked
    onto multiple lanes; the gene strip grows to fit. Returns
    ``(fig, coverage_axes)``; the gene axes (if any) is ``fig.axes[1]``.
    """
    if positions is None:
        positions = list(range(len(values)))
    if len(positions) != len(values):
        raise ValueError(
            f"positions and values must have the same length "
            f"(got {len(positions)} and {len(values)})"
        )

    if genes:
        lanes = _assign_lanes(genes)
        n_lanes = max(lanes) + 1
        gene_track_inches = max(0.5, 0.4 * n_lanes)
        gap_inches = 0.55
        avg_axes_h = (figsize[1] + gene_track_inches) / 2
        fig, (ax, gene_ax) = plt.subplots(
            2,
            1,
            figsize=(figsize[0], figsize[1] + gap_inches + gene_track_inches),
            dpi=dpi,
            gridspec_kw={
                "height_ratios": [figsize[1], gene_track_inches],
                "hspace": gap_inches / avg_axes_h,
            },
            sharex=True,
        )
        ax.tick_params(labelbottom=True)
    else:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        gene_ax = None
        lanes = None

    ax.fill_between(positions, values, step="mid", color=color, alpha=0.85, label=label)
    ax.plot(positions, values, drawstyle="steps-mid", color=color, linewidth=0.8)

    ax.set_xlim(positions[0], positions[-1])
    ax.set_ylim(bottom=0)
    ax.margins(x=0)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if label:
        ax.legend(loc="upper right", frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if gene_ax is not None:
        assert genes is not None  # invariant: gene_ax exists only when genes were given
        add_gene_track(
            gene_ax,
            genes,
            forward_color=forward_color,
            reverse_color=reverse_color,
            lanes=lanes,
        )
        gene_ax.set_xlabel("")
        boundary_xs = sorted({x for g in genes for x in (g.start, g.end)})
        for x in boundary_xs:
            ax.axvline(
                x,
                linestyle=(0, (1, 2)),
                color=guide_color,
                linewidth=0.7,
                alpha=0.8,
                zorder=0,
            )
            gene_ax.axvline(
                x,
                linestyle=(0, (1, 2)),
                color=guide_color,
                linewidth=0.7,
                alpha=0.8,
                zorder=0,
            )

    return fig, ax


def add_gene_track(
    gene_ax: Axes,
    genes: Sequence[Gene],
    *,
    forward_color: str = "#2563eb",
    reverse_color: str = "#dc2626",
    label: bool = True,
    fontsize: float = 8.0,
    arrow: bool = True,
    lanes: Sequence[int] | None = None,
) -> None:
    """Draw genes (with optional intron/exon structure) on a dedicated axes.

    Forward-strand (+) and reverse-strand (-) genes get different colors *and*
    a directional arrow tip on the terminal exon. Overlapping genes are
    auto-stacked into lanes (top-to-bottom). The axes is stripped of ticks,
    labels, and spines so it acts as a clean annotation strip.
    """
    if lanes is None:
        lanes = _assign_lanes(genes)
    n_lanes = max(lanes, default=0) + 1

    gene_ax.set_ylim(n_lanes, 0)  # lane 0 at top
    gene_ax.set_yticks([])
    gene_ax.tick_params(
        axis="x", which="both", length=0, labelbottom=False, labeltop=False
    )
    for spine in gene_ax.spines.values():
        spine.set_visible(False)

    # Arrow-tip width: a small fraction of the visible x-range so the tip
    # reads as consistent across the figure, regardless of gene length.
    xmin, xmax = gene_ax.get_xlim()
    arrow_w_base = max(1.0, (xmax - xmin) * 0.012)

    box_h = 0.55
    for gene, lane in zip(genes, lanes):
        gene_color = forward_color if gene.strand == "+" else reverse_color
        center_y = lane + 0.5
        box_y0 = center_y - box_h / 2

        # Thin intron line spanning the full gene
        gene_ax.plot(
            [gene.start, gene.end],
            [center_y, center_y],
            color=gene_color,
            linewidth=1.1,
            solid_capstyle="butt",
            zorder=1,
        )

        # Strand chevrons along the intron line; boxes will cover the ones
        # that fall inside exons, so visually they only show in introns.
        if arrow:
            span = gene.end - gene.start
            n_chev = max(1, int(span / 35))
            glyph = "›" if gene.strand == "+" else "‹"
            for k in range(n_chev):
                cx = gene.start + (k + 0.5) * span / n_chev
                gene_ax.text(
                    cx,
                    center_y,
                    glyph,
                    ha="center",
                    va="center",
                    color=gene_color,
                    fontsize=fontsize,
                    zorder=1,
                )

        # Exon shapes: terminal exon (in transcription direction) gets a
        # pointy arrow tip; the rest are plain rectangles.
        exons = sorted(gene.exons) if gene.exons else [(gene.start, gene.end)]
        terminal_idx = len(exons) - 1 if gene.strand == "+" else 0
        arrow_side: Literal["left", "right"] = "right" if gene.strand == "+" else "left"

        for i, (ex_start, ex_end) in enumerate(exons):
            if i == terminal_idx:
                tip_w = min(arrow_w_base, (ex_end - ex_start) * 0.5)
                verts = _exon_polygon(
                    ex_start, ex_end, box_y0, box_h,
                    arrow_side=arrow_side, arrow_width=tip_w,
                )
            else:
                verts = _exon_polygon(ex_start, ex_end, box_y0, box_h)
            gene_ax.add_patch(
                Polygon(verts, closed=True, facecolor=gene_color, edgecolor="none", zorder=2)
            )

        if label:
            gene_ax.text(
                (gene.start + gene.end) / 2,
                lane + 0.12,
                gene.name,
                ha="center",
                va="center",
                fontsize=fontsize - 0.5,
                color=gene_color,
                zorder=3,
                clip_on=False,
            )


def interactive_coverage_track(
    values: Sequence[float],
    positions: Sequence[int] | None = None,
    *,
    genes: Sequence[Gene] | None = None,
    selection_color: str = "#3b7dd8",
    crosshair_color: str = "#555555",
    dpi: int = 100,
    figsize: tuple[float, float] = (9.5, 2.6),
    **kwargs,
) -> tuple[Figure, Axes]:
    """Like :func:`coverage_track`, plus IGV-style mouse interactions.

    Designed for the ``ipympl`` Jupyter backend — activate it with
    ``%matplotlib widget`` in a cell *before* calling this function.

    Defaults to a lower ``dpi`` and a slightly smaller ``figsize`` than the
    static :func:`coverage_track` so the widget fits inside a notebook column
    without horizontal overflow. Override either if you want a bigger canvas.

    Interactions
    ------------
    * **Click + drag** on a track to zoom into the selected x-range.
    * **Double-click** anywhere on the figure to reset the zoom.
    * **Hover** to show a thin crosshair and a position readout that spans
      both the coverage and (if present) the gene track.

    Without an interactive backend (e.g. plain inline / Agg), the callbacks
    are still attached but won't fire; the figure renders as a normal static
    plot.
    """
    fig, ax = coverage_track(
        values, positions, genes=genes, dpi=dpi, figsize=figsize, **kwargs
    )

    # ipympl-specific chrome trimming so the widget fits in a notebook column.
    # These attributes only exist on the widget backend; ignore otherwise.
    canvas = fig.canvas
    for attr in ("header_visible", "footer_visible"):
        if hasattr(canvas, attr):
            try:
                setattr(canvas, attr, False)
            except Exception:
                pass
    if hasattr(canvas, "layout"):
        try:
            canvas.layout.width = "100%"
            canvas.layout.max_width = f"{int(figsize[0] * dpi)}px"
        except Exception:
            pass
    if hasattr(canvas, "capture_scroll"):
        try:
            canvas.capture_scroll = True
        except Exception:
            pass
    tracked_axes = list(fig.axes)
    original_xlim = ax.get_xlim()

    crosshairs = [
        a.axvline(
            x=original_xlim[0],
            color=crosshair_color,
            linewidth=0.7,
            alpha=0.0,
            zorder=10,
        )
        for a in tracked_axes
    ]
    cursor_text = ax.text(
        0.005,
        0.97,
        "",
        color=crosshair_color,
        fontsize=8,
        ha="left",
        va="top",
        transform=ax.transAxes,
        zorder=11,
    )

    def _on_move(event):
        if event.inaxes in tracked_axes and event.xdata is not None:
            x = event.xdata
            for line in crosshairs:
                line.set_xdata([x, x])
                line.set_alpha(0.7)
            cursor_text.set_text(f"x = {x:.0f}")
        else:
            for line in crosshairs:
                line.set_alpha(0.0)
            cursor_text.set_text("")
        fig.canvas.draw_idle()

    def _on_click(event):
        if event.dblclick and event.inaxes in tracked_axes:
            ax.set_xlim(original_xlim)
            fig.canvas.draw_idle()

    def _on_select(xmin, xmax):
        if xmax - xmin < 1:
            return
        ax.set_xlim(xmin, xmax)
        fig.canvas.draw_idle()

    span = SpanSelector(
        ax,
        _on_select,
        "horizontal",
        useblit=True,
        props={"alpha": 0.2, "facecolor": selection_color},
        interactive=False,
        minspan=1.0,
    )

    fig.canvas.mpl_connect("motion_notify_event", _on_move)
    fig.canvas.mpl_connect("button_press_event", _on_click)

    # Keep references alive so widgets aren't garbage-collected.
    fig._s2fview_widgets = (span,)  # type: ignore[attr-defined]

    return fig, ax


def _exon_polygon(
    start: float,
    end: float,
    y0: float,
    height: float,
    *,
    arrow_side: Literal["left", "right"] | None = None,
    arrow_width: float = 0.0,
) -> list[tuple[float, float]]:
    """Vertices for an exon: plain rectangle, or pentagon with a directional tip."""
    if not arrow_side or arrow_width <= 0:
        return [(start, y0), (end, y0), (end, y0 + height), (start, y0 + height)]
    arrow_width = min(arrow_width, end - start)
    if arrow_side == "right":
        body_end = end - arrow_width
        return [
            (start, y0),
            (body_end, y0),
            (end, y0 + height / 2),
            (body_end, y0 + height),
            (start, y0 + height),
        ]
    # left
    body_start = start + arrow_width
    return [
        (body_start, y0),
        (end, y0),
        (end, y0 + height),
        (body_start, y0 + height),
        (start, y0 + height / 2),
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
