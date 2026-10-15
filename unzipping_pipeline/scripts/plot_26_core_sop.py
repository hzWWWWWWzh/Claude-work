#!/usr/bin/env python3
"""Per-ROI panels and Luz histograms for the SOP-reprocessed 26-ROI cohort.

Panels follow SOP section 5.1.  Cycle 1 is red and cycle 2 is blue, and a Luz
histogram replaces the KDE peak signal as panel c.

House style: samples are small points with no connecting line, models are
lines, axes are four-sided with inward ticks, no background grid.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from reprocess_26_core_sop import (
    C_CYCLE1, C_CYCLE2, FULL_BP, HANDLE_BP, OUT,
    delta_nm_per_bp, x_ds,
)

PEAK_F_LO, PEAK_F_HI = 13.0, 32.7      # SOP section 6.1 peak window
HIST_BINS = np.arange(0.0, FULL_BP + 20.0, 20.0)
CYCLE_COLOUR = {1: C_CYCLE1, 2: C_CYCLE2}


def style_axes(ax):
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(direction="in", top=True, right=True, length=3.2, width=0.8,
                   labelsize=7)
    ax.grid(False)


def model_curve(force: np.ndarray, luz_bp: float, zero_f: float,
                full_f: float, full_bp: float = FULL_BP) -> np.ndarray:
    """Normalised extension of a state at constant Luz, in the 0-`full_bp` frame.

    `full_bp` is the construct length the batch is normalised to: 890 for the
    July/August constructs, 980 for YTrunk_2601_980bp (20261004).
    """
    ds = HANDLE_BP * (x_ds(force) - x_ds(np.full_like(force, zero_f)))
    e_full = (HANDLE_BP * (x_ds(np.array([full_f])) - x_ds(np.array([zero_f])))
              + full_bp * delta_nm_per_bp(np.array([full_f])))[0]
    model_scale = full_bp / e_full
    return (ds + luz_bp * delta_nm_per_bp(force)) * model_scale


def roi_page(pdf, roi: str, pts: pd.DataFrame, anc: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    a = anc[anc.roi_key == roi]
    cycles = sorted(pts.cycle.unique())

    # --- a) normalised FEC with the two dWLC state curves -----------------
    ax = axes[0]
    for cyc in cycles:
        col = CYCLE_COLOUR[cyc]
        for br, fill in (("up", True), ("down", False)):
            g = pts[(pts.cycle == cyc) & (pts.branch == br)]
            if g.empty:
                continue
            ax.plot(g.force_bin_pn, g.extension_normalized_nm, "o", ms=1.6,
                    mfc=col if fill else "none", mec=col, mew=0.5,
                    linestyle="none", alpha=0.85)
    row = a.iloc[0]
    grid = np.linspace(max(5.0, pts.force_bin_pn.min()), pts.force_bin_pn.max(), 250)
    for luz, ls in ((0.0, "--"), (FULL_BP, "-")):
        ax.plot(grid, model_curve(grid, luz, row.zero_force_pn, row.full_force_pn),
                ls, color="0.45", lw=0.8, zorder=0)
    ax.set_xlabel("Force (pN)", fontsize=8)
    ax.set_ylabel("Normalized extension (nm)", fontsize=8)
    ax.set_ylim(-250, 1150)
    style_axes(ax)
    ax.text(-0.22, 1.04, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")

    # --- b) Luz vs force ---------------------------------------------------
    ax = axes[1]
    for cyc in cycles:
        col = CYCLE_COLOUR[cyc]
        for br, fill in (("up", True), ("down", False)):
            g = pts[(pts.cycle == cyc) & (pts.branch == br)]
            g = g[g.force_bin_pn >= 8.0]
            if g.empty:
                continue
            ax.plot(g.force_bin_pn, g.luz_dwlc_bp, "o", ms=1.6,
                    mfc=col if fill else "none", mec=col, mew=0.5,
                    linestyle="none", alpha=0.85)
    for y in (0.0, 400.0, FULL_BP):
        ax.axhline(y, color="0.6", lw=0.6, ls=":", zorder=0)
    ax.set_xlabel("Force (pN)", fontsize=8)
    ax.set_ylabel("Luz (bp)", fontsize=8)
    ax.set_ylim(-150, 1100)
    style_axes(ax)
    ax.text(-0.22, 1.04, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")

    # --- c) Luz histogram, increase branch inside the SOP peak window ------
    ax = axes[2]
    any_counts = False
    for cyc in cycles:
        g = pts[(pts.cycle == cyc) & (pts.branch == "up")
                & pts.force_bin_pn.between(PEAK_F_LO, PEAK_F_HI)]
        v = g.luz_dwlc_bp.dropna()
        v = v[(v >= 0) & (v <= FULL_BP)]
        if v.empty:
            continue
        any_counts = True
        ax.hist(v, bins=HIST_BINS, histtype="step", lw=1.0,
                color=CYCLE_COLOUR[cyc], label=f"C{cyc}  n={len(v)}")
    if not any_counts:
        ax.text(0.5, 0.5, "no peak-eligible signal", ha="center", va="center",
                fontsize=7, color="0.4", transform=ax.transAxes)
    ax.set_xlabel("Luz (bp)", fontsize=8)
    ax.set_ylabel("Count", fontsize=8)
    ax.set_xlim(0, FULL_BP)
    if any_counts:
        ax.legend(fontsize=5.2, frameon=False, loc="upper right", handletextpad=0.4,
                  borderpad=0.2, labelspacing=0.2)
    style_axes(ax)
    ax.text(-0.22, 1.04, "c", transform=ax.transAxes, fontsize=9, fontweight="bold")

    bits = []
    for r in a.itertuples():
        note = ""
        if r.zero_shift_reason:
            note = f" (moved from {r.zero_nominal_pn:.0f})"
        bits.append(f"C{r.cycle}: 0 at {r.zero_force_pn:.1f} pN {r.zero_branch}{note}; "
                    f"890 at {r.full_force_pn:.0f} pN {r.full_branch}")
    fig.suptitle(f"{roi}   force x1.15 · Luz 0-890 bp · " + " | ".join(bits),
                 fontsize=6.4, y=1.02)
    fig.tight_layout()
    pdf.savefig(fig, dpi=300, bbox_inches="tight")
    plt.close(fig)


def cohort_histogram(pts: pd.DataFrame, path: Path):
    """Cohort Luz histogram, cycles kept separate per SOP section 6."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7))
    up = pts[(pts.branch == "up") & pts.force_bin_pn.between(PEAK_F_LO, PEAK_F_HI)]

    ax = axes[0]
    for cyc in (1, 2):
        v = up[up.cycle == cyc].luz_dwlc_bp.dropna()
        v = v[(v >= 0) & (v <= FULL_BP)]
        if v.empty:
            continue
        ax.hist(v, bins=HIST_BINS, histtype="step", lw=1.2,
                color=CYCLE_COLOUR[cyc],
                label=f"cycle {cyc}  n={len(v)}, {up[up.cycle==cyc].roi_key.nunique()} ROI")
    ax.set_xlabel("Luz (bp)", fontsize=8)
    ax.set_ylabel("Count (all ROI pooled)", fontsize=8)
    ax.set_xlim(0, FULL_BP)
    ax.legend(fontsize=6, frameon=False)
    style_axes(ax)
    ax.text(-0.16, 1.04, "a", transform=ax.transAxes, fontsize=9, fontweight="bold")

    # ROI-equal-weight version: each ROI contributes one normalised histogram,
    # so a densely sampled bead cannot dominate (SOP section 8).
    ax = axes[1]
    centres = 0.5 * (HIST_BINS[:-1] + HIST_BINS[1:])
    for cyc in (1, 2):
        rows = []
        for roi, g in up[up.cycle == cyc].groupby("roi_key"):
            v = g.luz_dwlc_bp.dropna()
            v = v[(v >= 0) & (v <= FULL_BP)]
            if len(v) < 5:
                continue
            h, _ = np.histogram(v, bins=HIST_BINS)
            rows.append(h / h.sum())
        if not rows:
            continue
        m = np.vstack(rows)
        ax.plot(centres, m.mean(axis=0), "-", lw=1.2, color=CYCLE_COLOUR[cyc],
                label=f"cycle {cyc}  {len(rows)} ROI")
        se = m.std(axis=0, ddof=1) / np.sqrt(len(rows)) if len(rows) > 1 else np.zeros_like(centres)
        ax.fill_between(centres, m.mean(axis=0) - se, m.mean(axis=0) + se,
                        color=CYCLE_COLOUR[cyc], alpha=0.18, lw=0)
    ax.set_xlabel("Luz (bp)", fontsize=8)
    ax.set_ylabel("Mean fraction per ROI", fontsize=8)
    ax.set_xlim(0, FULL_BP)
    ax.legend(fontsize=6, frameon=False)
    style_axes(ax)
    ax.text(-0.16, 1.04, "b", transform=ax.transAxes, fontsize=9, fontweight="bold")

    fig.suptitle(f"Luz histogram, increase branch {PEAK_F_LO}-{PEAK_F_HI} pN, "
                 f"20 bp bins · cycle 1 = red, cycle 2 = blue",
                 fontsize=7, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(path.with_suffix("." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    pts = pd.read_csv(OUT / "core26_normalized_points_sop.csv.gz", low_memory=False)
    anc = pd.read_csv(OUT / "per_cycle_anchors.csv")
    anc["zero_shift_reason"] = anc.zero_shift_reason.fillna("")

    rois = sorted(pts.roi_key.unique())
    pdf_path = OUT / "all_26_ROI_SOP_reprocessed_FEC.pdf"
    index = []
    with PdfPages(pdf_path) as pdf:
        for i, roi in enumerate(rois, start=1):
            roi_page(pdf, roi, pts[pts.roi_key == roi], anc)
            index.append({"pdf_page": i, "roi_key": roi,
                          "n_cycles": pts[pts.roi_key == roi].cycle.nunique()})
    pd.DataFrame(index).to_csv(OUT / "all_26_ROI_page_index.csv", index=False)
    print(f"{pdf_path.name}: {len(rois)} pages")

    cohort_histogram(pts, OUT / "cohort26_Luz_histogram")
    print("cohort26_Luz_histogram.{pdf,png,svg}")

    up = pts[(pts.branch == "up") & pts.force_bin_pn.between(PEAK_F_LO, PEAK_F_HI)]
    rows = []
    for (roi, cyc), g in up.groupby(["roi_key", "cycle"]):
        v = g.luz_dwlc_bp.dropna()
        v = v[(v >= 0) & (v <= FULL_BP)]
        if v.empty:
            continue
        h, _ = np.histogram(v, bins=HIST_BINS)
        rows.append({"roi_key": roi, "cycle": cyc, "n_points": len(v),
                     **{f"bin_{int(HIST_BINS[k])}_{int(HIST_BINS[k+1])}": int(h[k])
                        for k in range(len(h))}})
    pd.DataFrame(rows).to_csv(OUT / "per_roi_cycle_Luz_histogram.csv", index=False)
    print(f"per_roi_cycle_Luz_histogram.csv: {len(rows)} ROI-cycles")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
