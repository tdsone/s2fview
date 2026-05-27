from __future__ import annotations

import pathlib
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Literal

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from matplotlib.widgets import SpanSelector


def _register_bundled_fonts() -> None:
    """Register fonts shipped with the package so matplotlib can find them.

    We bundle Inter (variable) so the viewer looks the same on every machine
    without relying on the user having Inter pre-installed.
    """
    fonts_dir = pathlib.Path(__file__).parent / "fonts"
    if not fonts_dir.is_dir():
        return
    for font_path in fonts_dir.glob("*.ttf"):
        font_manager.fontManager.addfont(str(font_path))


_register_bundled_fonts()

# Prefer Inter, then a chain of macOS / common system fallbacks. Matplotlib
# walks the list and uses the first one it can resolve, silently falling
# back to the platform default otherwise.
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = [
    "Inter Variable",
    "Inter",
    "Avenir Next",
    "SF Pro Text",
    "Helvetica Neue",
    "Helvetica",
    *mpl.rcParams["font.sans-serif"],
]


# IGV-ish palette for DNA bases. Saturated enough to be unambiguous, soft
# enough that white letters on top read cleanly.
SEQUENCE_COLORS: dict[str, str] = {
    "A": "#3aa55a",  # green
    "C": "#3585c7",  # blue
    "G": "#e8a43c",  # orange
    "T": "#d65a4a",  # red
    "N": "#a8a8a8",  # gray fallback
}


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
    label: str | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str = "Coverage",
    figsize: tuple[float, float] = (10, 2.2),
    dpi: int = 200,
    guide_color: str = "#9a9a9a",
) -> tuple[Figure, Axes]:
    """Plot a coverage track, optionally with a DNA sequence and a gene track.

    When ``sequence`` and/or ``genes`` are provided, the figure becomes a
    stack of x-shared axes:

    * Coverage (top, with tick labels)
    * DNA sequence strip (colored cells per base; letters appear only when
      there's enough pixel room — they re-render lazily on zoom)
    * Gene annotation track (auto-stacked into lanes for overlapping genes)

    Returns ``(fig, coverage_axes)``; additional axes are available via
    ``fig.axes``.
    """
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

    has_sequence = sequence is not None
    has_genes = bool(genes)

    sizes_inches: list[float] = [figsize[1]]
    if has_sequence:
        sizes_inches.append(0.30)
    lanes: list[int] | None = None
    if has_genes:
        assert genes is not None
        lanes = _assign_lanes(genes)
        n_lanes = max(lanes) + 1
        sizes_inches.append(max(0.5, 0.4 * n_lanes))

    gene_ax: Axes | None = None
    seq_ax: Axes | None = None

    if len(sizes_inches) == 1:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        gap_inches = 0.45
        n_rows = len(sizes_inches)
        avg_axes_h = sum(sizes_inches) / n_rows
        fig = plt.figure(
            figsize=(figsize[0], sum(sizes_inches) + (n_rows - 1) * gap_inches),
            dpi=dpi,
        )
        gs = fig.add_gridspec(
            n_rows,
            1,
            height_ratios=sizes_inches,
            hspace=gap_inches / avg_axes_h,
        )
        ax = fig.add_subplot(gs[0])
        idx = 1
        if has_sequence:
            seq_ax = fig.add_subplot(gs[idx], sharex=ax)
            idx += 1
        if has_genes:
            gene_ax = fig.add_subplot(gs[idx], sharex=ax)
        ax.tick_params(labelbottom=True)

    ax.fill_between(positions, values, step="mid", color=color, alpha=0.85, label=label)
    ax.plot(positions, values, drawstyle="steps-mid", color=color, linewidth=0.8)

    # When a sequence is shown, widen xlim to match the cell extents so the
    # first/last bases aren't clipped in half.
    xlim_pad = 0.5 if has_sequence else 0
    ax.set_xlim(positions[0] - xlim_pad, positions[-1] + xlim_pad)
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

    if seq_ax is not None:
        assert sequence is not None
        add_sequence_track(seq_ax, sequence, start=int(positions[0]))

    if gene_ax is not None:
        assert genes is not None
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
            for guide_ax in (ax, gene_ax):
                guide_ax.axvline(
                    x,
                    linestyle=(0, (1, 2)),
                    color=guide_color,
                    linewidth=0.7,
                    alpha=0.8,
                    zorder=0,
                )

    return fig, ax


def add_sequence_track(
    seq_ax: Axes,
    sequence: str,
    start: int = 0,
    *,
    fontsize: float = 7.0,
    letter_color: str = "white",
) -> None:
    """Render a DNA sequence as a one-row colored strip on ``seq_ax``.

    The colored cells (one per base) are drawn as a single ``imshow`` so the
    cost is independent of sequence length. Individual base letters render
    on top *lazily* — they only appear when there's enough horizontal pixel
    room per base, and they're re-rendered for the visible window whenever
    the axes' xlim changes (e.g. on zoom).

    ``start`` is the genomic position of ``sequence[0]``.
    """
    n = len(sequence)
    rgb = _sequence_to_rgb_array(sequence)

    seq_ax.imshow(
        rgb,
        aspect="auto",
        extent=(start - 0.5, start + n - 0.5, 0, 1),
        interpolation="nearest",
        zorder=1,
    )
    seq_ax.set_ylim(0, 1)
    seq_ax.set_yticks([])
    seq_ax.tick_params(
        axis="x", which="both", length=0, labelbottom=False, labeltop=False
    )
    for spine in seq_ax.spines.values():
        spine.set_visible(False)

    # Match xlim to the sequence range if nothing has set it yet.
    if seq_ax.get_xlim() == (0.0, 1.0):
        seq_ax.set_xlim(start - 0.5, start + n - 0.5)

    letter_texts: list = []
    state = {"first_drawn": False}

    def _redraw_letters(*_):
        for txt in letter_texts:
            txt.remove()
        letter_texts.clear()

        bbox = seq_ax.get_window_extent()
        if bbox.width <= 0:
            return

        x0, x1 = seq_ax.get_xlim()
        pixels_per_base = bbox.width / max(1e-9, x1 - x0)
        # Need roughly the letter width plus a bit of padding.
        if pixels_per_base < fontsize * 1.2:
            return

        i_min = max(0, int(np.floor(x0 - start)))
        i_max = min(n, int(np.ceil(x1 - start)) + 1)

        for i in range(i_min, i_max):
            txt = seq_ax.text(
                start + i,
                0.5,
                sequence[i].upper(),
                ha="center",
                va="center",
                color=letter_color,
                fontsize=fontsize,
                fontweight="bold",
                zorder=2,
                clip_on=True,
            )
            letter_texts.append(txt)

    seq_ax.callbacks.connect("xlim_changed", _redraw_letters)

    def _on_first_draw(_event):
        if state["first_drawn"]:
            return
        state["first_drawn"] = True
        _redraw_letters()
        seq_ax.figure.canvas.draw_idle()

    seq_ax.figure.canvas.mpl_connect("draw_event", _on_first_draw)


def _sequence_to_rgb_array(sequence: str) -> np.ndarray:
    """Return a (1, len(sequence), 3) RGB array colored per nucleotide."""
    n = len(sequence)
    rgb = np.empty((1, n, 3), dtype=float)
    default = mcolors.to_rgb(SEQUENCE_COLORS["N"])
    palette = {b: mcolors.to_rgb(c) for b, c in SEQUENCE_COLORS.items()}
    for i, base in enumerate(sequence):
        rgb[0, i] = palette.get(base.upper(), default)
    return rgb


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
    show_toolbar: bool = True,
    **kwargs,
):
    """Like :func:`coverage_track`, plus IGV-style mouse interactions.

    Designed for the ``ipympl`` Jupyter backend — activate it with
    ``%matplotlib widget`` in a cell *before* calling this function.

    Defaults to a lower ``dpi`` and a slightly smaller ``figsize`` than the
    static :func:`coverage_track` so the widget fits inside a notebook column
    without horizontal overflow. Override either if you want a bigger canvas.

    Interactions
    ------------
    * **Toolbar** (when ``show_toolbar=True`` and ``ipywidgets`` is available):
      ``◀`` pan left, ``−`` zoom out, ``⌂`` reset, ``+`` zoom in, ``▶`` pan
      right.
    * **Click + drag** on a track to zoom into the selected x-range.
    * **Double-click** anywhere on the figure to reset the zoom.
    * **Hover** to show a thin crosshair and a position readout that spans
      both the coverage and (if present) the gene track.

    Return type
    -----------
    Under the widget backend with ``ipywidgets`` installed, returns an
    ``ipywidgets.VBox`` containing the toolbar + canvas — display it directly.
    Otherwise returns ``(fig, ax)`` so the static-backend path keeps working.
    """
    backend = mpl.get_backend().lower()
    is_widget = "ipympl" in backend or "widget" in backend or "nbagg" in backend

    # In widget mode, suppress matplotlib's auto-display so the figure only
    # appears inside our composed VBox.
    if is_widget and show_toolbar:
        plt.ioff()
    try:
        fig, ax = coverage_track(
            values, positions, genes=genes, dpi=dpi, figsize=figsize, **kwargs
        )
    finally:
        if is_widget and show_toolbar:
            plt.ion()

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

    # Throttle mouse-move updates to ~30 Hz. ipympl's blit() sends a diff
    # image whose region depends on canvas state and isn't always large
    # enough to clear the previous crosshair on the browser side (leaves
    # frozen vertical lines). Full draw_idle is correct but expensive; the
    # throttle keeps it cheap enough to feel responsive.
    import time

    _last_move = {"t": 0.0, "x": None}
    _move_period = 1.0 / 30.0

    def _on_move(event):
        now = time.monotonic()
        if now - _last_move["t"] < _move_period:
            return
        if event.inaxes in tracked_axes and event.xdata is not None:
            x = event.xdata
            if _last_move["x"] == x:
                return
            _last_move["x"] = x
            for line in crosshairs:
                line.set_xdata([x, x])
                line.set_alpha(0.7)
            cursor_text.set_text(f"x = {x:.0f}")
        else:
            if _last_move["x"] is None:
                return
            _last_move["x"] = None
            for line in crosshairs:
                line.set_alpha(0.0)
            cursor_text.set_text("")
        _last_move["t"] = now
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

    if is_widget and show_toolbar:
        toolbar = _make_zoom_toolbar(ax, fig, original_xlim)
        if toolbar is not None:
            try:
                import ipywidgets as widgets
            except ImportError:
                return fig, ax
            return widgets.VBox([toolbar, fig.canvas])

    return fig, ax


def _make_zoom_toolbar(ax: Axes, fig: Figure, original_xlim: tuple[float, float]):
    """Return an ipywidgets HBox of pan / zoom / reset buttons, or None."""
    try:
        import ipywidgets as widgets
    except ImportError:
        return None

    def _set_xlim(x0: float, x1: float) -> None:
        lo, hi = sorted(original_xlim)
        x0 = max(lo, x0)
        x1 = min(hi, x1)
        if x1 - x0 < 1:
            return
        ax.set_xlim(x0, x1)
        fig.canvas.draw_idle()

    def _zoom(factor: float) -> None:
        x0, x1 = ax.get_xlim()
        center = (x0 + x1) / 2
        new_half = (x1 - x0) / 2 * factor
        _set_xlim(center - new_half, center + new_half)

    def _pan(direction: float) -> None:
        x0, x1 = ax.get_xlim()
        shift = (x1 - x0) * 0.25 * direction
        # Clamp so we don't pan past the data edges.
        lo, hi = sorted(original_xlim)
        if shift < 0:
            shift = max(shift, lo - x0)
        else:
            shift = min(shift, hi - x1)
        _set_xlim(x0 + shift, x1 + shift)

    def _reset() -> None:
        ax.set_xlim(original_xlim)
        fig.canvas.draw_idle()

    btn_layout = widgets.Layout(width="38px", padding="0px")
    pan_l = widgets.Button(icon="caret-left", tooltip="Pan left", layout=btn_layout)
    zoom_out = widgets.Button(icon="search-minus", tooltip="Zoom out", layout=btn_layout)
    reset = widgets.Button(icon="home", tooltip="Reset zoom", layout=btn_layout)
    zoom_in = widgets.Button(icon="search-plus", tooltip="Zoom in", layout=btn_layout)
    pan_r = widgets.Button(icon="caret-right", tooltip="Pan right", layout=btn_layout)

    pan_l.on_click(lambda _b: _pan(-1))
    zoom_out.on_click(lambda _b: _zoom(2.0))
    reset.on_click(lambda _b: _reset())
    zoom_in.on_click(lambda _b: _zoom(0.5))
    pan_r.on_click(lambda _b: _pan(+1))

    return widgets.HBox(
        [pan_l, zoom_out, reset, zoom_in, pan_r],
        layout=widgets.Layout(justify_content="center"),
    )


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
