#!/usr/bin/env python3
"""Reusable WHZ publication-figure helpers.

The module intentionally distinguishes raw traces (markers only) from fits or
models (lines). It creates individual panels for later assembly in a vector
editor; it does not provide a manuscript-page compositor.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.figure import Figure
from matplotlib.lines import Line2D


SKILL_ROOT = Path(__file__).resolve().parents[1]
STYLE_PATH = SKILL_ROOT / "assets" / "default.mplstyle"

PALETTE = {
    # Categorical order; matches axes.prop_cycle in assets/default.mplstyle.
    # NPG Nature base (ggsci pal_npg), applied as the WHZ group palette.
    "red": "#E64B35",         # primary raw traces / force
    "blue": "#3C5488",        # second condition / cycle 2
    "green": "#00A087",       # teal green
    "orange": "#F39B7F",      # salmon / cycle 1 lighter accent
    "purple": "#AA3377",      # retained accent: nucleosome / boxplot points
    "grey": "#8491B4",        # neutral 6th category / background condition
    "cyan": "#4DBBD5",        # NPG cyan accent
    "bright_red": "#DC0000",  # NPG bright red accent
    "light_teal": "#91D1C2",  # NPG light teal accent
    "brown": "#7E6148",       # NPG brown accent
    "tan": "#B09C85",         # NPG tan accent
    "yellow": "#E6AB02",      # legacy alias from the 20260731 FEC deck; not part of the NPG base
    "pink": "#F39B7F",        # legacy alias for the salmon accent
    "light_grey": "#BFBFBF",  # reference lines, model families, controls
    "dark_grey": "#666666",   # neutral raw points (e.g. boxplot observations)
    "black": "#000000",       # fits, text, spines
}

# Final printed panel dimensions. Assemble exported panels in vector software.
PANEL_SIZES = {
    "full": (7.09, 4.25),   # 180 mm wide
    "half": (3.39, 2.65),   # 86 mm wide
    "square": (3.39, 3.39),
}


def apply_style() -> None:
    """Load the bundled WHZ rcParams into the current Matplotlib process."""

    plt.style.use(STYLE_PATH)


def new_panel(size: str = "half", **subplot_kwargs) -> tuple[Figure, Axes]:
    """Create one final-size panel, not a page-level multi-panel composite."""

    if size not in PANEL_SIZES:
        raise ValueError(f"Unknown size {size!r}; choose {sorted(PANEL_SIZES)}")
    apply_style()
    return plt.subplots(figsize=PANEL_SIZES[size], **subplot_kwargs)


def _paired_1d(x: Sequence[float], y: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    x_array = np.asarray(x)
    y_array = np.asarray(y)
    if x_array.ndim != 1 or y_array.ndim != 1:
        raise ValueError("x and y must both be one-dimensional")
    if x_array.size != y_array.size:
        raise ValueError("x and y must have equal length")
    return x_array, y_array


def style_axes(
    ax: Axes,
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
) -> Axes:
    """Apply the WHZ four-sided, inward-tick, grid-free data-axis style."""

    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(mpl.rcParams["axes.linewidth"])
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    return ax


def raw_trace(
    ax: Axes,
    x: Sequence[float],
    y: Sequence[float],
    *,
    color: str = PALETTE["blue"],
    label: str | None = None,
    marker: str = "o",
    marker_area: float = 6.0,
    filled: bool = False,
    alpha: float = 0.9,
    edge_width: float = 0.5,
    rasterized: bool = False,
    zorder: float = 2,
    **kwargs,
) -> PathCollection:
    """Plot raw trace samples as small circles or points without connecting lines."""

    x_array, y_array = _paired_1d(x, y)
    facecolors = color if filled or marker == "." else "none"
    artist = ax.scatter(
        x_array,
        y_array,
        s=marker_area,
        marker=marker,
        facecolors=facecolors,
        edgecolors=color,
        linewidths=edge_width,
        alpha=alpha,
        label=label,
        rasterized=rasterized,
        zorder=zorder,
        **kwargs,
    )
    artist.set_gid("whz-raw-trace")
    return artist


def fit_line(
    ax: Axes,
    x: Sequence[float],
    y: Sequence[float],
    *,
    color: str = PALETTE["black"],
    label: str | None = None,
    linewidth: float = 1.25,
    linestyle: str = "-",
    zorder: float = 3,
    **kwargs,
) -> Line2D:
    """Plot a fitted curve or primary model as a continuous line."""

    x_array, y_array = _paired_1d(x, y)
    (artist,) = ax.plot(
        x_array,
        y_array,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        label=label,
        zorder=zorder,
        **kwargs,
    )
    artist.set_gid("whz-fit-or-model")
    return artist


def model_family(
    ax: Axes,
    x: Sequence[float],
    curves: Iterable[Sequence[float]],
    *,
    color: str = PALETTE["light_grey"],
    linewidth: float = 0.7,
    linestyle: str = ":",
    alpha: float = 0.85,
) -> list[Line2D]:
    """Plot secondary/reference model curves without competing with the data."""

    artists = []
    for curve in curves:
        artists.append(
            fit_line(
                ax,
                x,
                curve,
                color=color,
                linewidth=linewidth,
                linestyle=linestyle,
                alpha=alpha,
                zorder=1,
            )
        )
    return artists


def summary_series(
    ax: Axes,
    x: Sequence[float],
    y: Sequence[float],
    *,
    yerr: Sequence[float] | None = None,
    color: str = PALETTE["blue"],
    label: str | None = None,
    connect: bool = True,
    marker: str = "o",
    jitter: float = 0.0,
    seed: int = 0,
    **kwargs,
):
    """Plot derived summaries; connect only when x has a meaningful order.

    Set ``jitter`` to a small x offset (in data units) when independent
    repeats would otherwise plot on top of each other at the same x.
    """

    x_array, y_array = _paired_1d(x, y)
    x_array = np.asarray(x_array, dtype=float)
    if jitter:
        rng = np.random.default_rng(seed)
        x_array = x_array + rng.uniform(-jitter, jitter, size=x_array.size)
    artist = ax.errorbar(
        x_array,
        y_array,
        yerr=None if yerr is None else np.asarray(yerr),
        color=color,
        marker=marker,
        markersize=3.2,
        markerfacecolor=color,
        markeredgewidth=0,
        linestyle="-" if connect else "none",
        linewidth=1.0,
        elinewidth=0.8,
        capsize=2,
        label=label,
        **kwargs,
    )
    return artist


def boxplot_with_points(
    ax: Axes,
    groups: Sequence[Sequence[float]],
    *,
    positions: Sequence[float] | None = None,
    colors: Sequence[str] | None = None,
    labels: Sequence[str] | None = None,
    seed: int = 0,
    jitter: float = 0.12,
) -> dict:
    """Draw transparent boxplots with every underlying observation visible.

    Positions must leave enough room for the box width: adjacent boxes must
    not overlap, otherwise the groups visually merge. Jitter spreads raw
    points inside each box so repeated cycles do not stack into a column.
    """

    arrays = [np.asarray(group) for group in groups]
    if not arrays or any(group.ndim != 1 for group in arrays):
        raise ValueError("groups must contain one-dimensional numeric sequences")
    if positions is None:
        positions = np.arange(1, len(arrays) + 1, dtype=float)
    positions = np.asarray(positions, dtype=float)
    if positions.size != len(arrays):
        raise ValueError("positions and groups must have equal length")
    if len(positions) > 1 and np.diff(np.sort(positions)).min() < 0.55:
        raise ValueError(
            "box positions are closer than the 0.55 box width; increase the "
            "spacing between category centres so boxes cannot overlap"
        )
    if colors is None:
        colors = [PALETTE["dark_grey"]] * len(arrays)
    if len(colors) != len(arrays):
        raise ValueError("colors and groups must have equal length")

    box = ax.boxplot(
        arrays,
        positions=positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        boxprops={"facecolor": "none", "edgecolor": "black", "linewidth": 0.9},
        whiskerprops={"color": "black", "linewidth": 0.8},
        capprops={"color": "black", "linewidth": 0.8},
        medianprops={"color": "black", "linewidth": 1.1},
    )
    rng = np.random.default_rng(seed)
    for position, values, color in zip(positions, arrays, colors):
        offsets = rng.uniform(-jitter, jitter, size=values.size)
        raw_trace(
            ax,
            np.full(values.size, position) + offsets,
            values,
            color=color,
            marker_area=7,
            alpha=0.75,
            zorder=2.5,
        )
    if labels is not None:
        if len(labels) != len(arrays):
            raise ValueError("labels and groups must have equal length")
        ax.set_xticks(positions, labels)
    return box


def reference_line(
    ax: Axes,
    *,
    x: float | None = None,
    y: float | None = None,
    color: str = PALETTE["light_grey"],
    linestyle: str = "--",
    linewidth: float = 0.8,
) -> Line2D:
    """Add one scientifically motivated reference line, never a full grid."""

    if (x is None) == (y is None):
        raise ValueError("provide exactly one of x or y")
    artist = (
        ax.axvline(x, color=color, linestyle=linestyle, linewidth=linewidth, zorder=0)
        if x is not None
        else ax.axhline(y, color=color, linestyle=linestyle, linewidth=linewidth, zorder=0)
    )
    artist.set_gid("whz-reference-line")
    return artist


def panel_label(
    ax: Axes,
    label: str,
    *,
    suffix: str = ")",
    x: float = -0.13,
    y: float = 1.02,
    fontsize: float = 10.5,
) -> None:
    """Add the default lowercase panel label, e.g. ``a)``."""

    artist = ax.text(
        x,
        y,
        f"{label}{suffix}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize,
        fontweight="bold",
        clip_on=False,
    )
    # Exempt from the in-axes text budget checked by validate_figure.
    artist.set_gid("whz-panel-label")


def use_common_limits(axes: Sequence[Axes], *, x: bool = True, y: bool = True) -> None:
    """Apply the union of current limits to panels showing comparable quantities."""

    axes = list(axes)
    if not axes:
        return
    if x:
        lower = min(axis.get_xlim()[0] for axis in axes)
        upper = max(axis.get_xlim()[1] for axis in axes)
        for axis in axes:
            axis.set_xlim(lower, upper)
    if y:
        lower = min(axis.get_ylim()[0] for axis in axes)
        upper = max(axis.get_ylim()[1] for axis in axes)
        for axis in axes:
            axis.set_ylim(lower, upper)


def validate_figure(fig: Figure, *, require_axis_labels: bool = True) -> list[str]:
    """Return actionable WHZ-style violations without changing the figure."""

    violations: list[str] = []
    for index, ax in enumerate(fig.axes, start=1):
        if ax.get_gid() == "whz-schematic":
            continue
        for name, spine in ax.spines.items():
            if not spine.get_visible():
                violations.append(f"axis {index}: {name} spine is hidden")
        visible_grid = [line for line in (*ax.get_xgridlines(), *ax.get_ygridlines()) if line.get_visible()]
        if visible_grid:
            violations.append(f"axis {index}: background grid lines are visible")
        if require_axis_labels:
            if not ax.get_xlabel().strip():
                violations.append(f"axis {index}: x-axis label with units is missing")
            if not ax.get_ylabel().strip():
                violations.append(f"axis {index}: y-axis label with units is missing")
        legend = ax.get_legend()
        if legend is not None and legend.get_frame().get_visible():
            violations.append(f"axis {index}: legend frame is visible")
        for line in ax.lines:
            if line.get_gid() == "whz-raw-trace" and line.get_linestyle() not in ("None", "none", "", " "):
                violations.append(f"axis {index}: raw trace has a connecting line")
        violations.extend(_text_budget_violations(ax, index))
        violations.extend(_mathtext_violations(ax, index))
    return violations


# Words belong in the caption. Anything beyond axis labels, tick labels, the
# panel letter and the legend is on a budget: at most MAX_TEXT_ELEMENTS blocks
# of at most MAX_TEXT_WORDS words each, and no complete sentences.
MAX_TEXT_ELEMENTS = 3
MAX_TEXT_WORDS = 6


def _mathtext_violations(ax: Axes, index: int) -> list[str]:
    """Mathtext is laid out glyph by glyph, so the label leaves the figure as one
    text object per character. Plain ASCII keeps it a single editable string."""

    problems: list[str] = []
    labelled = [("x label", ax.get_xlabel()), ("y label", ax.get_ylabel()),
                ("title", ax.get_title())]
    labelled += [("text", artist.get_text()) for artist in ax.texts]
    legend = ax.get_legend()
    if legend is not None:
        labelled += [("legend entry", entry.get_text()) for entry in legend.get_texts()]
    for kind, content in labelled:
        if content and content.count("$") >= 2:
            problems.append(
                f"axis {index}: {kind} uses mathtext ({content[:40]!r}) - write it in "
                "plain ASCII (G1, dG1, kBT, pN/nm) so it stays one text object")
    return problems


def _text_budget_violations(ax: Axes, index: int) -> list[str]:
    problems: list[str] = []
    budgeted = []
    for artist in ax.texts:
        content = artist.get_text().strip()
        if not content or artist.get_gid() == "whz-panel-label":
            continue
        budgeted.append(content)
        longest = max((len(line.split()) for line in content.splitlines()), default=0)
        if longest > MAX_TEXT_WORDS:
            problems.append(
                f"axis {index}: in-axes text is {longest} words, budget is "
                f"{MAX_TEXT_WORDS} - move it to the caption ({content.splitlines()[0][:48]!r})")
    if len(budgeted) > MAX_TEXT_ELEMENTS:
        problems.append(
            f"axis {index}: {len(budgeted)} free-standing text elements, budget is "
            f"{MAX_TEXT_ELEMENTS} - the panel is doing the caption's job")
    legend = ax.get_legend()
    # A legend identifies data series, so its length is set by the data, not by
    # taste; four is the practical ceiling before a reader stops matching them.
    if legend is not None and len(legend.get_texts()) > 4:
        problems.append(
            f"axis {index}: legend has {len(legend.get_texts())} entries, keep it to 4 "
            "or label the curves directly")
    return problems



_STYLE_ATTRIBUTE = re.compile(r'style="([^"]*)"')
_FONT_DECLARATION = re.compile(
    r"font:\s*((?:(?:italic|oblique|normal|bold|\d{3})\s+)*)([\d.]+)px\s+(.+)$"
)


def _split_font_shorthand(path: Path) -> int:
    """Rewrite CSS `font:` shorthand as separate font-family / font-size.

    Matplotlib emits `style="font: 9px 'Arial', 'Helvetica', ...'`. PowerPoint
    cannot parse that shorthand when the family list is quoted, drops the size
    with it and falls back to its own default, which renders every label about
    1.8x too large and clips the longest legend entry. Splitting the shorthand
    into explicit properties makes PowerPoint, Illustrator and browsers agree.

    The shorthand is found anywhere in the style attribute, not only at its
    start: a coloured label is written `style="fill: #e64b35; font: 8px ..."`,
    and an earlier version that anchored on `style="font:` left every coloured
    label at PowerPoint's default size. Style and weight may come in either
    order (`italic 700` or `700 italic`).
    """

    markup = path.read_text(errors="ignore")
    count = 0

    def rewrite(match: re.Match) -> str:
        nonlocal count
        declarations = [part.strip() for part in match.group(1).split(";") if part.strip()]
        out, changed = [], False
        for declaration in declarations:
            font = _FONT_DECLARATION.match(declaration)
            if not font:
                out.append(declaration)
                continue
            changed = True
            modifiers, size, family = font.groups()
            out += [f"font-family: {family.strip()}", f"font-size: {size}px"]
            for token in modifiers.split():
                if token in ("italic", "oblique"):
                    out.append(f"font-style: {token}")
                elif token != "normal":
                    out.append(f"font-weight: {token}")
        if not changed:
            return match.group(0)
        count += 1
        return 'style="' + "; ".join(out) + '"'

    rewritten = _STYLE_ATTRIBUTE.sub(rewrite, markup)
    if count:
        path.write_text(rewritten)
    return count


def _svg_text_is_editable(path: Path) -> list[str]:
    """Reject an SVG whose labels were written one x coordinate per glyph."""

    try:
        markup = path.read_text(errors="ignore")
    except OSError:
        return []
    problems = []
    for coords in re.findall(r'<tspan[^>]*\bx="([^"]*)"', markup):
        if len(coords.split()) > 1:
            problems.append(
                f"{path.name}: a tspan carries {len(coords.split())} x coordinates - "
                "the label will arrive as one text object per character; remove the "
                "mathtext that produced it")
            break
    return problems


def export_panel(
    fig: Figure,
    output_base: str | Path,
    *,
    formats: Sequence[str] = ("pdf", "svg", "png"),
    dpi: int = 600,
    validate: bool = True,
    require_axis_labels: bool = True,
    close: bool = True,
) -> list[Path]:
    """Export one panel in editable vector and high-resolution raster formats."""

    if validate:
        violations = validate_figure(fig, require_axis_labels=require_axis_labels)
        if violations:
            raise ValueError("WHZ figure validation failed:\n- " + "\n- ".join(violations))
    output_base = Path(output_base)
    if output_base.suffix:
        output_base = output_base.with_suffix("")
    output_base.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for file_format in formats:
        file_format = file_format.lower().lstrip(".")
        if file_format not in {"pdf", "svg", "png", "tif", "tiff"}:
            raise ValueError(f"unsupported output format: {file_format}")
        output = output_base.with_suffix(f".{file_format}")
        fig.savefig(output, format=file_format, dpi=dpi, facecolor="white")
        if file_format == "svg":
            _split_font_shorthand(output)
        if validate and file_format == "svg":
            svg_problems = _svg_text_is_editable(output)
            if svg_problems:
                raise ValueError("WHZ figure validation failed:\n- " + "\n- ".join(svg_problems))
        written.append(output)
    if close:
        plt.close(fig)
    return written
