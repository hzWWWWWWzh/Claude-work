#!/usr/bin/env python3
"""20261009: every ROI of every titration file, one page each, for the user to pick.

No pre-selection (user, 2026-10-11: "生成所有的FEC合集，我自己来挑选"). Every
ROI of the 17 titration files gets a page, in file and ROI order, with its
tracking verdict in the title; nothing is dropped.

Rules on these pages:
  * Zero = the MAXIMUM extension in the first cycle, the 10 pN step (user,
    2026-10-11: "命令10pN（第一个cycle）的最大extension为0点"): Hampel filter
    (+-0.5 s, 3 robust sigma) then the largest 1 s running median of that
    plateau, as in the 20261004 norm1028 rule for the maximum extension.
  * Drift: the bead goes back to the 5.00 mm rest (1.6 pN, fork closed) for
    15 s before every step, so the rest dwells sample the drift every 25 s.
    Extension = smoothed Z (median5 + mean5, as SOP 3.1) minus a line through
    the medians of the last 10 s of each rest. The SOP's line through the first
    and last second alone left 100-300 nm of mid-file drift (data_014 R039:
    a fully open fork that 'shrinks' from 1034 to 747 nm as the force rises),
    which is what produced the apparent stalls near 800 bp. A rest median more
    than 300 nm off its neighbours (tracking level) is skipped; with fewer than
    two good rests the SOP line is used and the title says so.
  * Luz is the dWLC model alone, no per-bead length scale:
    Luz = (ext - 1800 bp handle stretch from 10 pN) / (2 x_ss - x_ds), F >= 5 pN.
    Reference curves on the FEC: Luz 0 (closed), 314 / 460 / 511 / 657 (601
    site edges) and 980 bp (fully open).

    ~/.venvs/trace_editor_bm/bin/python fec_all_20261009.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter, uniform_filter1d

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_20261009 as p  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402
from process_tierA_sop import drift_corrected_nm  # noqa: E402
from reprocess_26_core_sop import (HANDLE_BP, MIN_DELTA_NM_PER_BP, MIN_INVERSION_FORCE,  # noqa: E402
                                   delta_nm_per_bp, x_ds)
from tracking_quality import file_tracking  # noqa: E402

OUT = p.DEST / "FEC_all"
ZERO_F_PN = 10.0
REST_TAIL_S = 10.0
REST_OUTLIER_NM = 300.0
SITE_EDGES = (314.0, 460.0, 511.0, 657.0)


def smoothed_nm(z_um: np.ndarray) -> np.ndarray:
    z = np.asarray(z_um, float)
    ok = np.isfinite(z)
    if ok.sum() < 2:
        return np.full_like(z, np.nan)
    if not ok.all():
        z = np.interp(np.arange(z.size), np.flatnonzero(ok), z[ok])
    return uniform_filter1d(median_filter(z, size=5, mode="nearest"), size=5, mode="nearest") * 1000.0


def rest_drift(run: dict, steps: list[dict], z_nm: np.ndarray) -> tuple[np.ndarray, str]:
    """Drift line through the closed-fork rest dwells (last 10 s of each)."""
    t, s = run["time_s"], run["shift_mm"]
    rest = np.abs(s - p.REST_MM) <= p.POS_TOL_MM
    pts = []
    for r in p.runs(rest):
        tail = r[t[r] >= t[r[-1]] - REST_TAIL_S]
        if tail.size >= 20 and t[r[-1]] - t[r[0]] >= 5.0:
            pts.append((float(np.median(t[tail])), float(np.nanmedian(z_nm[tail]))))
    if len(pts) >= 3:
        y = np.array([q[1] for q in pts])
        keep = np.ones(len(pts), bool)
        for i in range(len(pts)):
            nb = [y[j] for j in (i - 1, i + 1) if 0 <= j < len(pts)]
            keep[i] = abs(y[i] - np.median(nb)) <= REST_OUTLIER_NM
        pts = [q for q, k in zip(pts, keep) if k]
    if len(pts) < 2:
        return drift_corrected_nm(run["z_um_raw"], t), "SOP first/last-second line (< 2 good rests)"
    tx, ty = zip(*pts)
    return z_nm - np.interp(t, tx, ty), f"rest-dwell drift ({len(pts)} rests)"


def max_extension(t: np.ndarray, e: np.ndarray) -> float:
    """Hampel (+-0.5 s, 3 sigma), then the largest 1 s running median."""
    if e.size < 10 or not np.isfinite(e).any():
        return np.nan
    fps = (t.size - 1) / (t[-1] - t[0]) if t[-1] > t[0] else 40.0
    w = max(5, int(round(1.0 * fps)) | 1)
    s = pd.Series(e)
    med = s.rolling(w, center=True, min_periods=3).median()
    mad = (s - med).abs().rolling(w, center=True, min_periods=3).median() * 1.4826
    clean = s.where((s - med).abs() <= 3 * mad.replace(0, np.nan).fillna(np.inf))
    return float(clean.rolling(w, center=True, min_periods=max(3, w // 2)).median().max())


def model_ext(force: np.ndarray, luz: float) -> np.ndarray:
    return HANDLE_BP * (x_ds(force) - x_ds(np.full_like(force, ZERO_F_PN))) + luz * delta_nm_per_bp(force)


def model_luz(ext: np.ndarray, force: np.ndarray) -> np.ndarray:
    dx = delta_nm_per_bp(force)
    with np.errstate(divide="ignore", invalid="ignore"):
        luz = (ext - HANDLE_BP * (x_ds(force) - x_ds(np.full_like(force, ZERO_F_PN)))) / dx
    bad = (force < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP) | ~np.isfinite(luz)
    return np.where(bad, np.nan, luz)


def page(run, steps, roi, trk, colours, d):
    t = run["time_s"] - run["time_s"][0]
    f = run["force_pn"]
    z_nm = smoothed_nm(run["z_um"][roi])
    run["z_um_raw"] = run["z_um"][roi]
    ext, drift_note = rest_drift(run, steps, z_nm)
    pl1 = steps[0]["plateau"] if steps else np.array([], int)
    zero = max_extension(t[pl1], ext[pl1]) if pl1.size else np.nan
    e = ext - zero if np.isfinite(zero) else ext
    luz = model_luz(e, f) if np.isfinite(zero) else np.full_like(e, np.nan)
    lost = trk.get("status") == "lost" and np.isfinite(trk.get("t_loss_s", np.nan))
    t_loss = trk["t_loss_s"] - run["time_s"][0] if lost else np.inf

    step_of = np.zeros(t.size, int)
    for st in steps:
        for part in ("up", "plateau", "down"):
            step_of[st[part]] = st["step"]

    fig = plt.figure(figsize=(11.0, 5.4))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.55], height_ratios=[0.8, 1.2, 1.2],
                          hspace=0.10, wspace=0.22)
    axq = fig.add_subplot(gs[:, 0])
    axf = fig.add_subplot(gs[0, 1])
    axe = fig.add_subplot(gs[1, 1], sharex=axf)
    axl = fig.add_subplot(gs[2, 1], sharex=axf)

    after = t >= t_loss
    grey = (step_of == 0) | after
    axq.scatter(f[grey], e[grey], s=0.6, lw=0, color=PALETTE["light_grey"], rasterized=True)
    axe.scatter(t[grey], e[grey], s=0.6, lw=0, color=PALETTE["light_grey"], rasterized=True)
    for st in steps:
        m = (step_of == st["step"]) & ~after
        c = colours[round(st["force_pn"], 1)]
        axq.scatter(f[m], e[m], s=0.8, lw=0, color=c, rasterized=True)
        axe.scatter(t[m], e[m], s=0.8, lw=0, color=c, rasterized=True)
        axl.scatter(t[m], luz[m], s=0.8, lw=0, color=c, rasterized=True)
    if np.isfinite(zero):
        fg = np.linspace(MIN_INVERSION_FORCE, max(42.0, f.max()), 200)
        for L, ls in [(0.0, "--")] + [(x, ":") for x in SITE_EDGES] + [(p.CONSTRUCT_BP, "-")]:
            axq.plot(fg, model_ext(fg, L), ls, color=PALETTE["dark_grey"], lw=0.7, zorder=3)
    axf.plot(t, f, ".", ms=0.6, color=PALETTE["dark_grey"], rasterized=True)
    for lo, hi in p.NUCLEOSOME_SITES_BP:
        axl.axhspan(lo, hi, color=PALETTE["light_grey"], alpha=0.35, lw=0, zorder=0)
    for y, ls in ((0, "--"), (p.CONSTRUCT_BP, "-")):
        axl.axhline(y, color=PALETTE["light_grey"], ls=ls, lw=0.6, zorder=0)
    if lost:
        for ax in (axe, axl, axf):
            ax.axvline(t_loss, color=PALETTE["red"], lw=0.6, ls="--")

    lo, hi = np.nanpercentile(e[~after], [0.5, 99.5]) if (~after).any() else (-500, 1500)
    lo, hi = max(lo, -1500), min(hi, 2500)
    axq.set_ylim(min(lo, -200) - 50, max(hi, 900) + 50)
    axe.set_ylim(axq.get_ylim())
    axl.set_ylim(-300, 1.2 * p.CONSTRUCT_BP)
    style_axes(axq, xlabel="Force (pN)", ylabel="Extension from 10 pN zero (nm)")
    style_axes(axf, ylabel="Force (pN)")
    style_axes(axe, ylabel="Extension (nm)")
    style_axes(axl, xlabel="Time (s)", ylabel="Luz, dWLC (bp)")
    for ax in (axf, axe):
        ax.tick_params(labelbottom=False)
    status = trk.get("status", "unknown")
    extra = (f" · lost at {t_loss:.0f} s" if lost else "") + \
            (f" · level jump {trk.get('level_jump_um', 0):+.1f} um" if abs(trk.get("level_jump_um", 0)) > 0.5 else "")
    fig.suptitle(f"20261009 data_{d:03d} R{roi}   tracking {status}{extra} · zero = max ext at 10 pN "
                 f"(step 1) {zero:.0f} nm · {drift_note}", fontsize=7.5, y=0.995)
    return fig, {"data_id": d, "roi": f"R{roi}", "tracking": status, "t_loss_s": round(t_loss, 1) if lost else np.nan,
                 "level_jump_um": trk.get("level_jump_um", 0.0), "zero_10pN_nm": round(zero, 1) if np.isfinite(zero) else np.nan,
                 "drift": drift_note}


def main() -> int:
    p.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    titr, _ = p.titration_files()
    levels = sorted({10.0, 18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9, 46.0})
    colours = p.level_colours(levels)
    index = []
    for d in titr:
        run = p.load(d)
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{d:03d}.h5")
        pdf_path = OUT / f"WHZ_20261009_data_{d:03d}_FEC_all_{len(run['rois'])}_ROIs.pdf"
        with PdfPages(pdf_path) as pdf:
            for k, roi in enumerate(run["rois"], start=1):
                fig, row = page(run, steps, roi, trk.get(roi, {}), colours, d)
                pdf.savefig(fig, dpi=200, bbox_inches="tight")
                plt.close(fig)
                row.update(pdf=pdf_path.name, page=k)
                index.append(row)
        print(f"data_{d:03d}: {len(run['rois'])} pages -> {pdf_path.name}", flush=True)
    pd.DataFrame(index).to_csv(OUT / "FEC_all_page_index.csv", index=False)
    print(f"{len(index)} pages in {len(titr)} PDFs; index FEC_all_page_index.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
