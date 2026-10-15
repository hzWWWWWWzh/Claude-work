#!/usr/bin/env python3
"""Redraw all 59 processed ROIs in the WHZ figure style.

Covers the five PDFs in `Unzipping nucleoome/Proecessed Traces/`: the 26 core
ROIs of 20260731 plus the 33 ROIs held in the four image-only PDFs.  Data comes
from the SOP reprocessing (force x1.15, per-cycle anchors with the zero pulled
below any step inside its window, Luz normalised to 0-890 bp).

WHZ contract applied throughout: raw samples are small circles never joined by
lines, models are lines, all four spines with inward ticks, no grid, Arial at
9 pt, NPG palette with cycle 1 red and cycle 2 blue, panel labels a) b) c).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from whz_plot import (
    PALETTE, apply_style, export_panel, fit_line, panel_label, raw_trace,
    style_axes, validate_figure,
)

from reprocess_26_core_sop import FULL_BP
from plot_26_core_sop import HIST_BINS, PEAK_F_HI, PEAK_F_LO, model_curve

ROOT = Path(__file__).resolve().parent
DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces")
CORE = ROOT / "20260824_core26_SOP_reprocessed"
FOLDER = ROOT / "20260824_folder59_SOP"

CYCLE_COLOUR = {1: PALETTE["red"], 2: PALETTE["blue"]}   # C1 red, C2 blue


def _force_col(pts) -> str:
    """Raw-sample tables carry `force_pn`, binned tables `force_bin_pn`.

    Display always uses raw samples (see unzip_raw_points); the binned column is
    still accepted so older outputs can be regenerated unchanged.
    """
    return "force_pn" if "force_pn" in pts.columns else "force_bin_pn"
REF_GREY = PALETTE["light_grey"]
PAGE_SIZE = (7.09, 2.30)      # 180 mm wide row of three panels


def load_all() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Both SOP runs, tagged by which source PDF the ROI came from."""
    core = pd.read_csv(CORE / "core26_normalized_points_sop.csv.gz", low_memory=False)
    core_anc = pd.read_csv(CORE / "per_cycle_anchors.csv")
    core["date"] = 20260731
    core_anc["date"] = 20260731
    core["source_pdf"] = "all_26_ROI_core_branch_aware_FEC.pdf"
    core_anc["source_pdf"] = "all_26_ROI_core_branch_aware_FEC.pdf"
    core["anchor_valid"] = True

    other = pd.read_csv(FOLDER / "folder59_normalized_points_sop.csv.gz", low_memory=False)
    other_anc = pd.read_csv(FOLDER / "per_cycle_anchors.csv")

    shared = [c for c in ("date", "roi_key", "cycle", "branch", "force_bin_pn",
                          "extension_normalized_nm", "extension_mad_normalized_nm",
                          "luz_dwlc_bp", "luz_uncertainty_bp", "anchor_valid",
                          "source_pdf") if c in core.columns and c in other.columns]
    pts = pd.concat([core[shared], other[shared]], ignore_index=True)
    anc_cols = [c for c in ("date", "roi_key", "cycle", "zero_force_pn", "zero_branch",
                            "zero_nominal_pn", "zero_shift_reason", "full_force_pn",
                            "full_branch", "raw_anchor_span_nm", "anchor_valid",
                            "source_pdf") if c in core_anc.columns and c in other_anc.columns]
    anc = pd.concat([core_anc[anc_cols], other_anc[anc_cols]], ignore_index=True)
    anc["zero_shift_reason"] = anc.zero_shift_reason.fillna("")
    return pts[pts.anchor_valid], anc[anc.anchor_valid]


def draw_fec(ax, pts: pd.DataFrame, anc: pd.DataFrame, full_bp: float = FULL_BP):
    fcol = _force_col(pts)
    for cyc in sorted(pts.cycle.unique()):
        col = CYCLE_COLOUR.get(cyc, PALETTE["dark_grey"])
        for br, filled in (("up", True), ("down", False)):
            g = pts[(pts.cycle == cyc) & (pts.branch == br)]
            if g.empty:
                continue
            raw_trace(ax, g[fcol], g.extension_normalized_nm, color=col,
                      filled=filled, marker_area=3.0, alpha=0.8, edge_width=0.4,
                      rasterized=True)
    row = anc.iloc[0]
    grid = np.linspace(max(5.0, pts[fcol].min()), pts[fcol].max(), 300)
    # The two dWLC state boundaries are reference models, not fits to this trace.
    for luz, ls in ((0.0, "--"), (full_bp, "-")):
        fit_line(ax, grid, model_curve(grid, luz, row.zero_force_pn, row.full_force_pn,
                                       full_bp),
                 color=REF_GREY, linestyle=ls, linewidth=0.9, zorder=1)
    style_axes(ax, xlabel="Force (pN)", ylabel="Normalized extension (nm)")
    ax.set_ylim(-0.28 * full_bp, 1.30 * full_bp)
    panel_label(ax, "a")


def draw_luz(ax, pts: pd.DataFrame, full_bp: float = FULL_BP, bands=()):
    fcol = _force_col(pts)
    for cyc in sorted(pts.cycle.unique()):
        col = CYCLE_COLOUR.get(cyc, PALETTE["dark_grey"])
        for br, filled in (("up", True), ("down", False)):
            g = pts[(pts.cycle == cyc) & (pts.branch == br) & (pts[fcol] >= 8.0)]
            if g.empty:
                continue
            raw_trace(ax, g[fcol], g.luz_dwlc_bp, color=col, filled=filled,
                      marker_area=3.0, alpha=0.8, edge_width=0.4, rasterized=True)
    for y in (0.0, full_bp):
        ax.axhline(y, color=REF_GREY, lw=0.8, ls=":", zorder=1)
    # sequence features (e.g. 601 positioning sites) as shaded Luz bands
    for lo, hi in bands:
        ax.axhspan(lo, hi, color=REF_GREY, alpha=0.25, lw=0, zorder=0)
    style_axes(ax, xlabel="Force (pN)", ylabel=r"$L_{uz}$ (bp)")
    ax.set_ylim(-150, 1.24 * full_bp)
    panel_label(ax, "b")


def draw_hist(ax, pts: pd.DataFrame, full_bp: float = FULL_BP, bands=()):
    fcol = _force_col(pts)
    drew = False
    for cyc in sorted(pts.cycle.unique()):
        g = pts[(pts.cycle == cyc) & (pts.branch == "up")
                & pts[fcol].between(PEAK_F_LO, PEAK_F_HI)]
        v = g.luz_dwlc_bp.dropna()
        v = v[(v >= 0) & (v <= full_bp)]
        if v.empty:
            continue
        drew = True
        ax.hist(v, bins=np.arange(0.0, full_bp + 20.0, 20.0), histtype="step", lw=1.1,
                color=CYCLE_COLOUR.get(cyc, PALETTE["dark_grey"]),
                label=f"cycle {cyc}")
    if not drew:
        ax.text(0.5, 0.5, "no peak-eligible signal", ha="center", va="center",
                fontsize=8, color=PALETTE["dark_grey"], transform=ax.transAxes)
    else:
        ax.legend(loc="upper right")
    style_axes(ax, xlabel=r"$L_{uz}$ (bp)", ylabel="Count")
    for lo, hi in bands:
        ax.axvspan(lo, hi, color=REF_GREY, alpha=0.25, lw=0, zorder=0)
    ax.set_xlim(0, full_bp)
    panel_label(ax, "c")


def roi_page(pdf, date: int, roi: str, pts: pd.DataFrame, anc: pd.DataFrame, src: str):
    fig, axes = plt.subplots(1, 3, figsize=PAGE_SIZE)
    draw_fec(axes[0], pts, anc)
    draw_luz(axes[1], pts)
    draw_hist(axes[2], pts)
    bits = []
    for r in anc.itertuples():
        note = f" (from {r.zero_nominal_pn:.0f})" if r.zero_shift_reason else ""
        bits.append(f"C{r.cycle} zero {r.zero_force_pn:.1f} pN {r.zero_branch}{note}")
    fig.suptitle(f"{date}_{roi}   force ×1.15 · $L_{{uz}}$ 0–890 bp · " + "; ".join(bits)
                 + f"   [{src}]", fontsize=7, y=1.03)
    fig.tight_layout(w_pad=1.6)
    pdf.savefig(fig, dpi=600, bbox_inches="tight")
    plt.close(fig)


def cohort_panels(pts: pd.DataFrame, dest: Path):
    """Cohort Luz histogram, exported as standalone panels per the WHZ workflow."""
    up = pts[(pts.branch == "up") & pts.force_bin_pn.between(PEAK_F_LO, PEAK_F_HI)]
    centres = 0.5 * (HIST_BINS[:-1] + HIST_BINS[1:])

    apply_style()
    fig, ax = plt.subplots(figsize=(3.39, 2.65))
    for cyc in (1, 2):
        v = up[up.cycle == cyc].luz_dwlc_bp.dropna()
        v = v[(v >= 0) & (v <= FULL_BP)]
        if v.empty:
            continue
        n_roi = up[up.cycle == cyc].groupby(["date", "roi_key"]).ngroups
        ax.hist(v, bins=HIST_BINS, histtype="step", lw=1.2, color=CYCLE_COLOUR[cyc],
                label=f"cycle {cyc} ({n_roi} ROI)")
    style_axes(ax, xlabel=r"$L_{uz}$ (bp)", ylabel="Count")
    ax.set_xlim(0, FULL_BP)
    ax.legend(loc="upper center")
    panel_label(ax, "a")
    export_panel(fig, dest / "WHZ_cohort_Luz_histogram_pooled", dpi=600)

    fig, ax = plt.subplots(figsize=(3.39, 2.65))
    for cyc in (1, 2):
        rows = []
        for _, g in up[up.cycle == cyc].groupby(["date", "roi_key"]):
            v = g.luz_dwlc_bp.dropna()
            v = v[(v >= 0) & (v <= FULL_BP)]
            if len(v) < 5:
                continue
            h, _ = np.histogram(v, bins=HIST_BINS)
            rows.append(h / h.sum())
        if not rows:
            continue
        m = np.vstack(rows)
        mean = m.mean(axis=0)
        se = m.std(axis=0, ddof=1) / np.sqrt(len(rows)) if len(rows) > 1 else np.zeros_like(mean)
        ax.fill_between(centres, mean - se, mean + se, color=CYCLE_COLOUR[cyc],
                        alpha=0.20, lw=0, zorder=1)
        fit_line(ax, centres, mean, color=CYCLE_COLOUR[cyc], linewidth=1.2,
                 label=f"cycle {cyc} ({len(rows)} ROI)")
    style_axes(ax, xlabel=r"$L_{uz}$ (bp)", ylabel="Mean fraction per ROI")
    ax.set_xlim(0, FULL_BP)
    ax.legend(loc="upper center")
    panel_label(ax, "b")
    export_panel(fig, dest / "WHZ_cohort_Luz_histogram_ROI_weighted", dpi=600)


def main() -> int:
    apply_style()
    DEST.mkdir(parents=True, exist_ok=True)
    pts, anc = load_all()
    keys = anc[["date", "roi_key", "source_pdf"]].drop_duplicates().sort_values(
        ["source_pdf", "date", "roi_key"])
    print(f"ROIs to draw: {len(keys)}")

    out_pdf = DEST / "WHZ_59_ROI_FEC_Luz_histogram.pdf"
    index = []
    with PdfPages(out_pdf) as pdf:
        for r in keys.itertuples():
            p = pts[(pts.date == r.date) & (pts.roi_key == r.roi_key)]
            a = anc[(anc.date == r.date) & (anc.roi_key == r.roi_key)]
            if p.empty or a.empty:
                continue
            roi_page(pdf, int(r.date), r.roi_key, p, a, r.source_pdf)
            index.append({"pdf_page": len(index) + 1, "date": int(r.date),
                          "roi_key": r.roi_key, "source_pdf": r.source_pdf,
                          "n_cycles": int(p.cycle.nunique())})
    pd.DataFrame(index).to_csv(DEST / "WHZ_59_ROI_page_index.csv", index=False)
    print(f"{out_pdf.name}: {len(index)} pages")

    cohort_panels(pts, DEST)
    print("cohort panels exported (pdf/svg/png @600 dpi)")

    # deterministic style check on a representative page
    fig, axes = plt.subplots(1, 3, figsize=PAGE_SIZE)
    r = keys.iloc[0]
    draw_fec(axes[0], pts[(pts.date == r.date) & (pts.roi_key == r.roi_key)],
             anc[(anc.date == r.date) & (anc.roi_key == r.roi_key)])
    draw_luz(axes[1], pts[(pts.date == r.date) & (pts.roi_key == r.roi_key)])
    draw_hist(axes[2], pts[(pts.date == r.date) & (pts.roi_key == r.roi_key)])
    problems = validate_figure(fig)
    plt.close(fig)
    print("validate_figure:", problems if problems else "clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
