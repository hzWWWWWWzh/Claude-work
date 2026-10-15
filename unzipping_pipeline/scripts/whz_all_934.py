#!/usr/bin/env python3
"""All 934 ROIs of 20260731, SOP-processed and drawn in the WHZ figure style.

One PDF per data file, matching the 15 existing pptx decks.  Every ROI gets a
page, including the ones that cannot be normalised: those show the raw
drift-corrected extension instead of a fabricated 0-890 coordinate, and say so
on the page, because SOP section 4.1 excludes a cycle with span < 200 nm from
normalised FEC and Luz rather than rescaling it anyway.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from whz_plot import PALETTE, apply_style, panel_label, raw_trace, style_axes

from process_tierA_sop import (
    FORCE_SCALE, bin_branch, drift_corrected_nm, force_from_shift, split_branches,
)
from reprocess_26_core_sop import (
    ANCHOR_EXCEPTIONS, ANCHOR_HALF_WINDOW, DEFAULT_ANCHOR, FULL_BP, HANDLE_BP,
    MIN_DELTA_NM_PER_BP, MIN_INVERSION_FORCE, delta_nm_per_bp, estimate_anchor,
    first_step_force, x_ds,
)
from plot_26_core_sop import HIST_BINS, PEAK_F_HI, PEAK_F_LO, model_curve
from plot_whz_59 import (
    CYCLE_COLOUR, PAGE_SIZE, REF_GREY, _force_col, draw_fec, draw_hist, draw_luz,
)

RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/ unzipping nucleosome/20260731")
DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/"
            "Proecessed Traces/WHZ_20260731_all_ROI")
DATA_FILES = [3, 5, 9, 11, 13, 14, 16, 17, 18, 19, 21, 23, 25, 26, 28]
MIN_SPAN_NM = 200.0     # SOP 4.1 anchor_valid floor
# Anchor windows on raw samples. The zero sits just below the first unzipping
# step, which can cut the closed state off at any force, so its window stays
# tight. The full anchor sits on the down branch, which is one fully open state
# from Fmax all the way down to ~15 pN and varies smoothly with force, so a wide
# window extrapolated by Theil-Sen costs precision, not bias. Fast 47 s ramps
# carry only ~1 sample per pN near 38 pN, which is what these limits allow for.
RAW_ZERO_MIN_POINTS, RAW_ZERO_MAX_REACH_PN = 5, 2.5
RAW_FULL_MIN_POINTS, RAW_FULL_MAX_REACH_PN = 3, 4.0


def process_file(data_id: int, raw_dir: Path | None = None,
                 bin_pn: float | None = None,
                 calibration: str = "sop_20260724_x1.15") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw HDF5 -> binned, anchored, normalised points for every ROI in one file."""
    path = (raw_dir or RAW) / f"data_{data_id:03d}.h5"
    with h5py.File(path, "r") as h5:
        traces = h5["Traces"]
        time_s = np.asarray(traces["Time (s)"], dtype=float)
        shift_mm = np.asarray(traces["Shift (mm)"], dtype=float)
        rois = sorted(traces["Z (um)"].keys())
        z_all = {r: np.asarray(traces["Z (um)"][r], dtype=float) for r in rois}
    force_all = force_from_shift(shift_mm, calibration)
    branches = split_branches(time_s, shift_mm)

    points, anchors = [], []
    for roi in rois:
        roi_key = f"D{data_id:03d}_R{roi}"
        ext = drift_corrected_nm(z_all[roi], time_s)
        binned = []
        for cycle, direction, seg in branches:
            t = (bin_branch(force_all[seg], ext[seg], bin_pn) if bin_pn
                 else bin_branch(force_all[seg], ext[seg]))
            if t.empty:
                continue
            t.insert(0, "branch", direction)
            t.insert(0, "cycle", cycle)
            binned.append(t)
        g_all = (pd.concat(binned, ignore_index=True) if binned else
                 pd.DataFrame(columns=["cycle", "branch", "force_bin_pn",
                                       "force_median_pn", "extension_raw_nm",
                                       "extension_mad_raw_nm", "n_points"]))
        g_all.insert(0, "roi_key", roi_key)
        # Anchors are estimated on RAW samples (already median5+mean5 smoothed by
        # drift_corrected_nm), not on bins. On fast ramps the >=2-points-per-bin
        # rule leaves no bins near 14 or 38 pN at all, so a binned anchor was
        # really measured several pN away (99% of cycles on 20261004). The
        # frame reuses the binned column names so the estimators are unchanged.
        g_raw = pd.concat([pd.DataFrame({"cycle": cycle, "branch": direction,
                                         "force_bin_pn": force_all[seg],
                                         "extension_raw_nm": ext[seg]})
                           for cycle, direction, seg in branches], ignore_index=True)

        for cyc in sorted(g_raw.cycle.unique()):
            g = g_all[g_all.cycle == cyc]
            r = g_raw[g_raw.cycle == cyc]
            zf, zb, ff, fb = ANCHOR_EXCEPTIONS.get(roi_key, DEFAULT_ANCHOR)
            zbr, fbr = r[r.branch == zb], r[r.branch == fb]
            z_fb = f_fb = ""
            if zbr.empty or zbr.force_bin_pn.min() > zf + ANCHOR_HALF_WINDOW:
                alt = "up" if zb == "down" else "down"
                if not r[r.branch == alt].empty:
                    z_fb, zb, zbr = f"{zb}->{alt}", alt, r[r.branch == alt]
            if fbr.empty or fbr.force_bin_pn.max() < ff - ANCHOR_HALF_WINDOW:
                alt = "up" if fb == "down" else "down"
                cand = r[r.branch == alt]
                if not cand.empty and cand.force_bin_pn.max() >= ff - ANCHOR_HALF_WINDOW:
                    f_fb, fb, fbr = f"{fb}->{alt}", alt, cand

            rec = dict(data_id=data_id, roi_key=roi_key, cycle=int(cyc),
                       anchor_valid=False, reason="", zero_force_pn=np.nan,
                       zero_branch="", zero_nominal_pn=np.nan, zero_shift_reason="",
                       full_force_pn=np.nan, full_branch="", raw_anchor_span_nm=np.nan,
                       zero_branch_fallback=z_fb, full_branch_fallback=f_fb)
            p = g.copy()
            p["extension_normalized_nm"] = np.nan
            p["luz_dwlc_bp"] = np.nan
            p["luz_uncertainty_bp"] = np.nan

            if zbr.empty or fbr.empty:
                rec["reason"] = "anchor branch missing"
            else:
                noise = (float(np.nanmedian(g.extension_mad_raw_nm))
                         if len(g) and np.isfinite(g.extension_mad_raw_nm).any() else 10.0)
                step_f = first_step_force(zbr, noise, zf, lag=5)
                # Raw samples: ask for 5 points and allow up to 2.5 pN of reach.
                # On a 253 s ramp the +-0.6 pN window already holds dozens, so
                # nothing changes; on a 47 s ramp 38 pN carries ~1 sample/pN.
                za = estimate_anchor(zbr, zf, avoid_above_pn=step_f,
                                     min_points=RAW_ZERO_MIN_POINTS,
                                     max_reach_pn=RAW_ZERO_MAX_REACH_PN)
                ff_eff = (min(ff, float(fbr.force_bin_pn.max()) - ANCHOR_HALF_WINDOW)
                          if fbr.force_bin_pn.max() < ff + ANCHOR_HALF_WINDOW else ff)
                fa = estimate_anchor(fbr, ff_eff, min_points=RAW_FULL_MIN_POINTS,
                                     max_reach_pn=RAW_FULL_MAX_REACH_PN)
                if za is None or fa is None:
                    rec["reason"] = ("zero anchor: too few samples" if za is None
                                     else "full anchor: too few samples")
                else:
                    span = fa.value_nm - za.value_nm
                    rec.update(zero_force_pn=za.force_pn, zero_branch=za.branch,
                               zero_nominal_pn=za.shifted_from_pn,
                               zero_shift_reason=za.shift_reason,
                               full_force_pn=fa.force_pn, full_branch=fa.branch,
                               raw_anchor_span_nm=round(span, 2),
                               # absolute anchor heights, so a later step can
                               # renormalise with a span taken from elsewhere
                               z_zero_raw_nm=round(za.value_nm, 2),
                               z_full_raw_nm=round(fa.value_nm, 2))
                    if span < MIN_SPAN_NM:
                        rec["reason"] = f"span {span:.0f} nm < {MIN_SPAN_NM:.0f} nm"
                    else:
                        rec["anchor_valid"] = True
                        scale = FULL_BP / span
                        p["extension_normalized_nm"] = (p.extension_raw_nm - za.value_nm) * scale
                        f = p.force_bin_pn.values
                        dx = delta_nm_per_bp(f)
                        ds = HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, za.force_pn)))
                        e_full = (HANDLE_BP * (x_ds(np.array([fa.force_pn]))
                                               - x_ds(np.array([za.force_pn])))
                                  + FULL_BP * delta_nm_per_bp(np.array([fa.force_pn])))[0]
                        ms = FULL_BP / e_full if e_full else np.nan
                        with np.errstate(divide="ignore", invalid="ignore"):
                            luz = (p.extension_normalized_nm.values / ms - ds) / dx
                            unc = np.abs(p.extension_mad_raw_nm.values * scale / ms / dx)
                        bad = ((f < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP)
                               | ~np.isfinite(luz))
                        p["luz_dwlc_bp"] = np.where(bad, np.nan, luz)
                        p["luz_uncertainty_bp"] = np.where(bad, np.nan, unc)
            p["anchor_valid"] = rec["anchor_valid"]
            anchors.append(rec)
            points.append(p)
    return (pd.concat(points, ignore_index=True) if points else pd.DataFrame(),
            pd.DataFrame(anchors))


def raw_page(fig, axes, pts: pd.DataFrame, anc: pd.DataFrame):
    """Fallback page: raw drift-corrected extension, nothing normalised."""
    ax = axes[0]
    fcol = _force_col(pts)
    for cyc in sorted(pts.cycle.unique()):
        col = CYCLE_COLOUR.get(cyc, PALETTE["dark_grey"])
        for br, filled in (("up", True), ("down", False)):
            g = pts[(pts.cycle == cyc) & (pts.branch == br)]
            if g.empty:
                continue
            raw_trace(ax, g[fcol], g.extension_raw_nm, color=col, filled=filled,
                      marker_area=3.0, alpha=0.8, edge_width=0.4, rasterized=True)
    style_axes(ax, xlabel="Force (pN)", ylabel="Drift-corrected extension (nm)")
    panel_label(ax, "a")
    for ax, letter in ((axes[1], "b"), (axes[2], "c")):
        ax.text(0.5, 0.5, "not normalizable\n(see page title)", ha="center", va="center",
                fontsize=8, color=PALETTE["dark_grey"], transform=ax.transAxes)
        style_axes(ax, xlabel="Force (pN)" if letter == "b" else r"$L_{uz}$ (bp)",
                   ylabel=r"$L_{uz}$ (bp)" if letter == "b" else "Count")
        panel_label(ax, letter)


def main() -> int:
    apply_style()
    DEST.mkdir(parents=True, exist_ok=True)
    summary, all_anchors = [], []
    for data_id in DATA_FILES:
        pts, anc = process_file(data_id)
        if pts.empty:
            print(f"data_{data_id:03d}: no data")
            continue
        rois = sorted(pts.roi_key.unique())
        out = DEST / f"WHZ_20260731_data_{data_id:03d}_{len(rois)}_ROIs.pdf"
        index = []
        with PdfPages(out) as pdf:
            for roi in rois:
                p = pts[pts.roi_key == roi]
                a = anc[anc.roi_key == roi]
                fig, axes = plt.subplots(1, 3, figsize=PAGE_SIZE)
                valid = a[a.anchor_valid]
                if len(valid):
                    pv = p[p.anchor_valid]
                    draw_fec(axes[0], pv, valid)
                    draw_luz(axes[1], pv)
                    draw_hist(axes[2], pv)
                    bits = [f"C{r.cycle} zero {r.zero_force_pn:.1f} pN {r.zero_branch}"
                            + (f" (from {r.zero_nominal_pn:.0f})" if r.zero_shift_reason else "")
                            for r in valid.itertuples()]
                    dropped = a[~a.anchor_valid]
                    if len(dropped):
                        bits += [f"C{r.cycle} excluded: {r.reason}" for r in dropped.itertuples()]
                    note = "; ".join(bits)
                else:
                    raw_page(fig, axes, p, a)
                    note = "; ".join(f"C{r.cycle} excluded: {r.reason}" for r in a.itertuples())
                fig.suptitle(f"20260731_{roi}   force ×{FORCE_SCALE} · "
                             f"$L_{{uz}}$ 0–890 bp · {note}", fontsize=7, y=1.03)
                fig.tight_layout(w_pad=1.6)
                pdf.savefig(fig, dpi=600, bbox_inches="tight")
                plt.close(fig)
                index.append({"pdf_page": len(index) + 1, "roi_key": roi,
                              "n_valid_cycles": int(a.anchor_valid.sum()),
                              "n_cycles": len(a)})
        pd.DataFrame(index).to_csv(
            DEST / f"WHZ_20260731_data_{data_id:03d}_page_index.csv", index=False)
        anc["data_id"] = data_id
        all_anchors.append(anc)
        nv = int(anc.groupby("roi_key").anchor_valid.any().sum())
        summary.append({"data_id": data_id, "n_roi": len(rois),
                        "roi_with_valid_cycle": nv,
                        "n_cycles": len(anc), "valid_cycles": int(anc.anchor_valid.sum()),
                        "pdf": out.name})
        print(f"data_{data_id:03d}: {len(rois)} ROI -> {out.name} "
              f"({nv} normalizable, {int(anc.anchor_valid.sum())}/{len(anc)} cycles)",
              flush=True)

    s = pd.DataFrame(summary)
    s.to_csv(DEST / "WHZ_20260731_summary.csv", index=False)
    pd.concat(all_anchors, ignore_index=True).to_csv(DEST / "per_cycle_anchors.csv", index=False)
    print(f"\nTOTAL {s.n_roi.sum()} ROI, {s.roi_with_valid_cycle.sum()} normalizable, "
          f"{s.valid_cycles.sum()}/{s.n_cycles.sum()} cycles")
    (DEST / "manifest.json").write_text(json.dumps({
        "sop": "20260811_unzipnucleosomedata_complete_processing_pipeline.docx",
        "date": 20260731, "data_files": DATA_FILES,
        "n_roi": int(s.n_roi.sum()), "force_scale": FORCE_SCALE,
        "full_length_bp": FULL_BP, "min_anchor_span_nm": MIN_SPAN_NM,
        "style": "whz-figure-style; cycle 1 red #E64B35, cycle 2 blue #3C5488",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
