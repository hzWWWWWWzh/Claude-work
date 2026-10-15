#!/usr/bin/env python3
"""20261009 ladder: the user's 96 selected tethers, merged; raw, FEC, normalised trace and Luz peak map.

Selection: the pages the user kept in the FEC_all PDFs (2026-10-11), read back
from the backup copy _user_selection_backup_96/ (selection_96_tethers.csv).

Rules (user, 2026-10-11):
  * 0 point = the maximum extension in the first cycle, the 10 pN step
    ("命令10pN（第一个cycle）的最大extension为0点"): fec_all_20261009.max_extension
    (Hampel +-0.5 s, 3 robust sigma, then the largest 1 s running median of
    that plateau). Drift: line through the 1.6 pN rests (fec_all.rest_drift).
  * 1028 point = the maximum extension (same estimator) on the plateau of the
    highest-force step the bead SURVIVED ("如果最后一次这个珠子已经断了（飞了）就不要
    用它的最终长度作为归一化最长extension，用上一个力"). A step counts only if
      - the bead is not lost before the end of the rest that follows it;
      - its plateau never runs (1 s running median) more than FLY_NM above or
        BELOW_NM (1000 nm) below the zero (fly-away / tracking failure, data_029
        R048); 300-1000 nm below is only flagged (data_046 R051: its 10 pN step
        sits 250 nm above the 18-30 pN plateaus, another LUT level?);
      - the 1.6 pN rest after it is back within REST_OK_NM of the closed
        baseline (a bead that flew does not come back); a rise > END_JUMP_NM in
        the plateau's last second counts as flying off only when there is no
        rest after it to show the bead came back (data_042 R016 vs data_003
        R023, whose +902 nm jump at 30 pN was the fork crossing the sites).
    The first step that fails ends the ladder: it and every later step are
    skipped and the step before is used.
  * normalised = 1028 * (extension - zero) / (end - zero); 601 sites at
    314-460 / 511-657.
  * FEC kept (user: "保留fec，但是增加luz的peak map"): extension vs force with the
    dWLC curves for Luz 0 / 314 / 460 / 511 / 657 / 980 bp from the 10 pN zero.
  * Luz, model B (user, 2026-10-11: "用模型B重新算Luz和峰图"): each opened bp adds
    2 x_ss(F) - the closed trunk dangles beyond the fork and is not in the force
    path (the Y geometry of the user's schematic), unlike the earlier model A
    (2 x_ss - x_ds, ~1.45x more bp per nm). Luz_B = (ext - zero - 1800 bp
    handle stretch from 10 pN) / (2 x_ss(F)), ABSOLUTE bp from the 10 pN zero,
    not rescaled per tether (rescaling the 1028 point to 1028 would make the
    model choice matter by only ~3 %). MODEL=A reproduces the earlier pages.
  * Luz peak map (basic-unzipping-analysis rules) on Luz_B:
    histogram of the samples at F >= 13 pN of the usable steps (before any
    loss), 5 bp bins, square-root normalised counts, adaptive local-background
    peak calling (160 bp window, 20 % quantile, min distance 80 bp, <= 5 peaks),
    the skill's make_histogram.

Outputs (WHZ_20261009/):
  WHZ_20261009_selected_96_FEC_pages.pdf   the user's kept FEC pages, merged as they were
  WHZ_20261009_selected_96_norm1028_modelB.pdf   one page per tether, bookmarked
  norm1028_selected_96_modelB.csv                 zero, reference step, end (nm and Luz_B),
                                                  span, skipped steps, Luz_B peaks, flags
  luz_modelB/fig_luzB_peakmap_pooled.(pdf|png)    the 96 tethers pooled (each tether weight 1),
                                                  all forces and by force group

    .venv/bin/python ladder_norm1028_20261009.py
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path.home() / ".claude/skills/basic-unzipping-analysis/scripts"))
import process_20261009 as p  # noqa: E402
import fec_all_20261009 as F  # noqa: E402
import basic_unzipping_analysis as BUA  # noqa: E402

import fitz  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

NORM_BP = 1028.0
FLY_NM = 2000.0        # a 980 bp fork fully open at 40 pN is ~1.2 um above the zero
BELOW_NM = 1000.0      # a plateau this far below the zero: tracking failure (ends the ladder)
BELOW_FLAG_NM = 300.0  # below the zero by this much: flagged only (zero or LUT level?)
BRIEF_MAX_S = 1.5      # a maximum held shorter than this is flagged for the user
END_JUMP_NM = 600.0    # a rise this big in the last second of a plateau: the bead flies
REST_OK_NM = 300.0     # the rest after a survived step is back at the closed baseline
NORM_LIM = (-300.0, 1300.0)
LUZ_LIM = (-200.0, 1250.0)
HIST_F_MIN = 13.0
SEL_CSV = p.DEST / "selection_96_tethers.csv"
BACKUP = p.DEST / "_user_selection_backup_96"
SITE_EDGES = (314.0, 460.0, 511.0, 657.0)
import os as _os
from reprocess_26_core_sop import HANDLE_BP, MIN_INVERSION_FORCE, x_ds, x_ss  # noqa: E402
MODEL = _os.environ.get("MODEL", "B")
# SCALE (user, 2026-10-11: "使用我的方法吧"): the peak map is made on the user's own
# per-tether 0-1028 normalisation and sits next to that trace; SCALE=luz gives the
# model-Luz layout of the previous version
SCALE = _os.environ.get("SCALE", "norm1028")
SUFFIX = "_peakmap" if SCALE == "norm1028" else ("_modelB" if MODEL == "B" else "")


def gain_nm_per_bp(f):
    f = np.asarray(f, float)
    return 2 * x_ss(f) if MODEL == "B" else 2 * x_ss(f) - x_ds(f)


def model_ext(f, luz):
    """Extension above the 10 pN zero of a fork opened by luz bp (FEC curves)."""
    f = np.asarray(f, float)
    return HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, F.ZERO_F_PN))) + luz * gain_nm_per_bp(f)


def model_luz(e_rel, f):
    f = np.asarray(f, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        luz = (e_rel - HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, F.ZERO_F_PN)))) / gain_nm_per_bp(f)
    return np.where((f < MIN_INVERSION_FORCE) | ~np.isfinite(luz), np.nan, luz)


def peak_args() -> argparse.Namespace:
    """The skill's default (adaptive) histogram and peak settings, 1028 full length."""
    return argparse.Namespace(
        full_unzip_bp=NORM_BP, hist_force_min_pn=HIST_F_MIN, hist_bin_bp=5.0, hist_smooth_bins=9,
        peak_method="adaptive", peak_local_window_bp=160.0, peak_baseline_quantile=0.20,
        peak_noise_sigma=1.5, peak_min_height=0.05, peak_min_prominence=0.04, peak_min_luz_bp=100.0,
        peak_min_distance_bp=80.0, max_peaks=5, energy_guide_window_bp=40.0)


def max_and_time(t: np.ndarray, e: np.ndarray) -> tuple[float, float]:
    v = F.max_extension(t, e)
    w = max(5, int(round((t.size - 1) / (t[-1] - t[0]))) | 1) if t[-1] > t[0] else 5
    rm = pd.Series(e).rolling(w, center=True, min_periods=3).median().values
    return v, float(t[int(np.nanargmax(rm))]) if np.isfinite(rm).any() else np.nan


def bead_scale(run, steps, roi, trk):
    t, s_ = run["time_s"], run["shift_mm"]
    z_nm = F.smoothed_nm(run["z_um"][roi])
    run["z_um_raw"] = run["z_um"][roi]
    ext, drift_note = F.rest_drift(run, steps, z_nm)
    lost = trk.get("status") == "lost" and np.isfinite(trk.get("t_loss_s", np.nan))
    t_loss = trk["t_loss_s"] if lost else np.inf
    rec = {"zero_nm": np.nan, "t_zero": np.nan, "ref_step": np.nan, "ref_force_pn": np.nan,
           "end_nm": np.nan, "t_end": np.nan, "span_nm": np.nan, "drift": drift_note,
           "t_loss": t_loss, "skipped": "", "usable_steps": [], "flags": []}
    pl1 = steps[0]["plateau"]
    if t[pl1[-1]] >= t_loss:
        rec["flags"].append("lost before the end of the 10 pN step: no zero")
        return ext, rec
    zero, tz = max_and_time(t[pl1], ext[pl1])
    rec.update(zero_nm=zero, t_zero=tz)
    rec["usable_steps"].append(steps[0]["step"])
    ok_steps, skipped, below = [], [], []
    for k, st in enumerate(steps):
        if st["step"] == 1:
            continue
        pl = st["plateau"]
        t_after = t[steps[k + 1]["up"][0]] if k + 1 < len(steps) and steps[k + 1]["up"].size else t[-1]
        why = ""
        if t_after >= t_loss or t[pl[-1]] >= t_loss:
            why = f"lost at {t_loss - t[0]:.0f} s"
        else:
            tp, e = t[pl], ext[pl] - zero
            fps = (pl.size - 1) / (tp[-1] - tp[0]) if tp[-1] > tp[0] else 40.0
            # sustained levels only (1 s running median): single spikes do not count
            r = pd.Series(e).rolling(max(3, int(round(1.0 * fps))), center=True, min_periods=3).median().values
            last, before = tp >= tp[-1] - 1.0, tp < tp[-1] - 2.0
            rest = (np.abs(s_ - p.REST_MM) <= p.POS_TOL_MM) & (t > tp[-1]) & (t <= t_after)
            rest_med = np.nan
            if rest.sum() >= 20:
                rest &= t >= t[rest].max() - 10.0
                rest_med = float(np.nanmedian(ext[rest]))
            end_jump = (np.nanmax(r[last]) - np.nanmedian(r[before])) if before.sum() >= 5 else 0.0
            if np.nanmax(r) > FLY_NM:
                why = f"runs > {FLY_NM:.0f} nm above the zero (flies)"
            elif np.nanmin(r) < -BELOW_NM:
                why = f"falls {-np.nanmin(r):.0f} nm below the zero for > 1 s (tracking failure)"
            elif np.isfinite(rest_med) and abs(rest_med) > REST_OK_NM:
                why = f"rest after it at {rest_med:+.0f} nm, not back to closed (bead flew or tether changed)"
            elif not np.isfinite(rest_med) and end_jump > END_JUMP_NM:
                # a big jump at the very end is a fork advance if the bead comes
                # back closed afterwards (data_003 R023); with no rest after it
                # to show that, it is taken as the bead flying off (data_042 R016)
                why = f"jumps +{end_jump:.0f} nm at the end with no rest after it (flies off)"
        if why:
            skipped.append(f"{st['force_pn']:.1f} pN ({why})")
            break
        if np.nanmin(r) < -BELOW_FLAG_NM:
            below.append(f"{st['force_pn']:.1f} pN {np.nanmin(r):.0f} nm")
        ok_steps.append(st)
        rec["usable_steps"].append(st["step"])
    rec["skipped"] = "; ".join(skipped)
    if not ok_steps:
        rec["flags"].append("no step survived after 10 pN")
        return ext, rec
    ref = max(ok_steps, key=lambda s: s["force_pn"])
    pl = ref["plateau"]
    end, te = max_and_time(t[pl], ext[pl])
    rec.update(ref_step=ref["step"], ref_force_pn=round(ref["force_pn"], 1), end_nm=end, t_end=te,
               span_nm=end - zero)
    if skipped:
        rec["flags"].append("1028 from the step before: skipped " + rec["skipped"])
    if below:
        rec["flags"].append("plateau(s) below the 10 pN zero (zero on another LUT level?): " + ", ".join(below))
    # how long the fork held the maximum level (within 100 nm of it)
    held = float(np.sum((ext[pl] >= end - 100) & np.isfinite(ext[pl])) * np.median(np.diff(t[pl])))
    rec["max_held_s"] = held
    if held < BRIEF_MAX_S:
        rec["flags"].append(f"1028 point held only {held:.1f} s (a brief excursion - check)")
    if not np.isfinite(rec["span_nm"]) or rec["span_nm"] <= 0:
        rec["flags"].append("span <= 0: not normalised")
    elif rec["span_nm"] < 150:
        rec["flags"].append(f"small span {rec['span_nm']:.0f} nm (x{NORM_BP / rec['span_nm']:.1f})")
    elif rec["span_nm"] > 1450:
        rec["flags"].append(f"span {rec['span_nm']:.0f} nm > a full 980 bp unzip: LUT level or double tether?")
    if abs(trk.get("level_jump_um", 0.0)) > 0.5:
        rec["flags"].append(f"LUT level jump {trk['level_jump_um']:+.1f} um somewhere in the file")
    return ext, rec


def luz_abs(ext, force, rec):
    """Absolute Luz (MODEL B by default) from the 10 pN zero, not rescaled."""
    return model_luz(ext - rec["zero_nm"], force)


def draw(run, steps, roi, ext, rec, colours, title):
    t = run["time_s"] - run["time_s"][0]
    f = run["force_pn"]
    t_loss = rec["t_loss"] - run["time_s"][0]
    after = t >= t_loss
    step_of = np.zeros(t.size, int)
    for st in steps:
        for part in ("up", "plateau", "down"):
            step_of[st[part]] = st["step"]
    usable = np.isin(step_of, rec["usable_steps"]) & ~after
    # v2 (user, 2026-10-11: the peak map next to the 0-1028 trace "很wired"): the
    # Luz peak map now sits next to a Luz trace on the SAME axis (model B bp);
    # the user's 0-1028 trace is kept above it
    if SCALE == "norm1028":
        # raw + FEC; the user's 0-1028 trace + its peak map on the same axis; force
        fig = plt.figure(figsize=(7.09, 6.2))
        gs = fig.add_gridspec(3, 2, height_ratios=[1.5, 1.6, 0.7], width_ratios=[3.0, 1.55],
                              hspace=0.22, wspace=0.12)
        axr = fig.add_subplot(gs[0, 0])
        axn = fig.add_subplot(gs[1, 0], sharex=axr)
        axf = fig.add_subplot(gs[2, 0], sharex=axr)
        axq = fig.add_subplot(gs[0, 1])
        axp = fig.add_subplot(gs[1, 1], sharey=axn)
        axt = fig.add_subplot(gs[2, 1])
        axl = None
    else:
        fig = plt.figure(figsize=(7.09, 7.4))
        gs = fig.add_gridspec(4, 2, height_ratios=[1.4, 1.3, 1.4, 0.6], width_ratios=[3.0, 1.55],
                              hspace=0.22, wspace=0.12)
        axr = fig.add_subplot(gs[0, 0])
        axn = fig.add_subplot(gs[1, 0], sharex=axr)
        axl = fig.add_subplot(gs[2, 0], sharex=axr)
        axf = fig.add_subplot(gs[3, 0], sharex=axr)
        axq = fig.add_subplot(gs[0, 1])
        axt = fig.add_subplot(gs[1, 1])
        axp = fig.add_subplot(gs[2, 1], sharey=axl)
        axs = fig.add_subplot(gs[3, 1])
        axs.axis("off")
    luz_all = luz_abs(ext, f, rec) if np.isfinite(rec["zero_nm"]) else np.full_like(ext, np.nan)
    span = rec["span_nm"] if np.isfinite(rec["span_nm"]) and rec["span_nm"] > 0 else np.nan
    nrm = NORM_BP * (ext - rec["zero_nm"]) / span
    grey = (step_of == 0) | after
    for ax, y in ((axr, ext), (axn, nrm)):
        ax.scatter(t[grey], y[grey], s=0.5, lw=0, color=PALETTE["light_grey"], rasterized=True)
    axq.scatter(f[grey], (ext - rec["zero_nm"])[grey], s=0.5, lw=0, color=PALETTE["light_grey"],
                rasterized=True)
    for st in steps:
        m = (step_of == st["step"]) & ~after
        c = colours[round(st["force_pn"], 1)] if st["step"] in rec["usable_steps"] else PALETTE["grey"]
        axr.scatter(t[m], ext[m], s=0.7, lw=0, color=c, rasterized=True)
        axn.scatter(t[m], nrm[m], s=0.7, lw=0, color=c, rasterized=True)
        if axl is not None:
            axl.scatter(t[m], luz_all[m], s=0.7, lw=0, color=c, rasterized=True)
        axq.scatter(f[m], (ext - rec["zero_nm"])[m], s=0.7, lw=0, color=c, rasterized=True)
    for lo, hi in p.NUCLEOSOME_SITES_BP:
        for ax in (axn, axl, axp):
            if ax is not None:
                ax.axhspan(lo, hi, color=PALETTE["light_grey"], alpha=0.45, lw=0, zorder=0)
    for y in (0.0, NORM_BP):
        axn.axhline(y, color="black", lw=0.6, ls="--", zorder=4)
        if SCALE == "norm1028":
            axp.axhline(y, color="black", lw=0.6, ls="--", zorder=4)
    if axl is not None:
        for y in (0.0, 971.0):                   # Luz axis: closed and the trunk end
            for ax in (axl, axp):
                ax.axhline(y, color="black", lw=0.6, ls="--", zorder=4)
    if np.isfinite(rec["zero_nm"]):
        axr.axhline(rec["zero_nm"], color="black", lw=0.6, ls="--", zorder=4)
        axr.plot(rec["t_zero"] - run["time_s"][0], rec["zero_nm"], "^", ms=5, mfc="black",
                 mec="white", zorder=6)
        fg = np.linspace(5.0, 42.0, 200)
        for L, ls in [(0.0, "--")] + [(x, ":") for x in SITE_EDGES] + [(p.CONSTRUCT_BP, "-")]:
            axq.plot(fg, model_ext(fg, L), ls, color=PALETTE["dark_grey"], lw=0.6, zorder=3)
    if np.isfinite(rec["end_nm"]):
        axr.axhline(rec["end_nm"], color="black", lw=0.6, ls="--", zorder=4)
        axr.plot(rec["t_end"] - run["time_s"][0], rec["end_nm"], "o", ms=5, mfc="black",
                 mec="white", zorder=6)
        axq.plot(rec["ref_force_pn"], rec["end_nm"] - rec["zero_nm"], "o", ms=5, mfc="black",
                 mec="white", zorder=6)
    # Luz peak map
    peaks = []
    if np.isfinite(span):
        luz = nrm if SCALE == "norm1028" else luz_abs(ext, f, rec)
        m = usable & np.isfinite(luz)
        rec["luz_samples"] = (luz[m], f[m])
        centers, sqrt_norm, smooth, base, sig, peaks, _ = BUA.make_histogram(luz[m], f[m], peak_args())
        peaks = [x for x in peaks if abs(x - NORM_BP) > 1e-6]      # the full marker is drawn anyway
        axp.barh(centers, sqrt_norm, height=5.0, color=PALETTE["light_grey"], lw=0)
        axp.plot(base, centers, color=PALETTE["dark_grey"], lw=0.7, ls="--")
        axp.plot(sig, centers, color=PALETTE["blue"], lw=0.8)
        for x in peaks:
            axp.axhline(x, xmin=0.75, xmax=1.0, color=PALETTE["red"], lw=1.0)
            axp.text(1.02, x, f"{x:.0f}", transform=axp.get_yaxis_transform(), color=PALETTE["red"],
                     fontsize=6, va="center", ha="left")
        axp.set_xlim(0, 1.05)
    axf.plot(t, f, ".", ms=0.5, color=PALETTE["dark_grey"], rasterized=True)
    if np.isfinite(t_loss):
        for ax in (axr, axn, axl, axf):
            if ax is not None:
                ax.axvline(t_loss, color=PALETTE["red"], lw=0.7, ls="--")
    ok = ~after & np.isfinite(ext) & (step_of > 0)
    if ok.any():
        lo, hi = np.nanpercentile(ext[ok], [1, 99])
        lo = min(lo, rec["zero_nm"] if np.isfinite(rec["zero_nm"]) else lo) - 100
        hi = max(hi, rec["end_nm"] if np.isfinite(rec["end_nm"]) else hi) + 100
        axr.set_ylim(lo, max(hi, lo + 1300))
        axq.set_ylim(lo - rec["zero_nm"], max(hi, lo + 1300) - rec["zero_nm"])
    axn.set_ylim(*NORM_LIM)
    if axl is not None:
        axl.set_ylim(*LUZ_LIM)
    axq.set_xlim(0, 45)
    style_axes(axr, ylabel="Extension (nm)")
    style_axes(axn, ylabel=f"Normalised (bp)\n1028 = {span:.0f} nm" if np.isfinite(span) else "Normalised (bp)")
    if axl is not None:
        style_axes(axl, ylabel=f"Luz, model {MODEL} (bp)")
    style_axes(axf, xlabel="Time (s)", ylabel="Force (pN)")
    style_axes(axq, xlabel="Force (pN)", ylabel="Ext. from zero (nm)")
    style_axes(axp, xlabel="sqrt counts (norm.)",
               ylabel="Normalised (0-1028)" if SCALE == "norm1028" else f"Luz, model {MODEL} (bp)")
    axq.yaxis.set_label_position("right"); axq.yaxis.tick_right()
    axp.yaxis.set_label_position("right"); axp.yaxis.tick_right()
    for ax in (axr, axn, axl):
        if ax is not None:
            ax.tick_params(labelbottom=False)
    axt.axis("off")
    if SCALE != "norm1028" and np.isfinite(rec.get("end_luz", np.nan)) and rec["end_luz"] < 0.9 * 971 and not any(
            "not a full unzip" in w for w in rec["flags"]):
        rec["flags"].append(f"1028 point = Luz_{MODEL} {rec['end_luz']:.0f} bp: not a full unzip; the 0-1028 "
                            f"scale stretches this tether x{NORM_BP / max(rec['end_luz'], 1):.2f}")
    lines = [f"0: max ext at 10 pN (step 1) {rec['zero_nm']:.0f} nm",
             (f"1028: max ext at {rec['ref_force_pn']:.1f} pN {rec['end_nm']:.0f} nm"
              if np.isfinite(rec["ref_step"]) else "1028: none"),
             f"span {rec['span_nm']:.0f} nm; {rec['drift']}",
             (f"1028 point = Luz_{MODEL} {rec['end_luz']:.0f} bp" if np.isfinite(rec.get("end_luz", np.nan)) else ""),
             ("Peaks, 0-1028 scale (F >= 13 pN): " if SCALE == "norm1028" else f"Luz_{MODEL} peaks (F >= 13 pN): ")
             + (", ".join(f"{x:.0f}" for x in peaks) if peaks else "none")]
    lines = [w for ln in lines for w in textwrap.wrap(ln, 42, subsequent_indent="  ")]
    axt.text(0.0, 0.72, "\n".join(lines), transform=axt.transAxes, va="top", ha="left", fontsize=6.3)
    if rec["flags"]:
        fig.text(0.01, -0.02, "\n".join(textwrap.fill("! " + w, 150) for w in rec["flags"]),
                 va="top", ha="left", fontsize=6.5, color=PALETTE["red"])
    fig.suptitle(f"{title}   (raw + FEC; " + ("normalised 0-1028 + its peak map on the same axis; "
                 if SCALE == "norm1028" else f"normalised 0-1028; Luz model {MODEL} + its peak map; ")
                 + "grey = after a loss / skipped step)",
                 fontsize=7.5, y=0.995)
    return fig, peaks


def pooled_peakmap(pooled):
    """Luz peak map of all selected tethers together (each tether weight 1), all
    forces and by force group, with the skill's adaptive local-background peaks."""
    out = p.DEST / ("norm1028_by_force" if SCALE == "norm1028" else f"luz_model{MODEL}")
    out.mkdir(exist_ok=True)
    a = peak_args()
    edges = np.arange(0.0, 1100.0 + 5.0, 5.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    groups = [("all forces >= 13 pN", 13.0, 99.0), ("18-27.9 pN", 17.0, 28.5),
              ("30-34.6 pN", 29.0, 35.5), ("37.1-46 pN", 36.0, 99.0)]
    fig, axes = plt.subplots(len(groups), 1, figsize=(7.09, 8.0), sharex=True)
    rows = []
    for ax, (lab, lo, hi) in zip(axes, groups):
        h = np.zeros(centers.size)
        n_t = 0
        for name, luz, f in pooled:
            m = (f >= max(lo, HIST_F_MIN)) & (f < hi) & np.isfinite(luz)
            if m.sum() < 20:
                continue
            c, _ = np.histogram(luz[m], bins=edges)
            h += c / m.sum()                       # each tether weight 1
            n_t += 1
        sq = np.sqrt(h)
        smax = sq.max() if sq.max() > 0 else 1.0
        sm = BUA.centered_mean(sq, a.hist_smooth_bins)
        base, sig = BUA.local_peak_signal(centers, sm, a.peak_local_window_bp, a.peak_baseline_quantile)
        peaks = BUA.find_histogram_peaks(centers, sm, sig, sm / smax, "adaptive", a.peak_min_height,
                                         a.peak_min_prominence * smax, a.peak_noise_sigma,
                                         a.peak_min_distance_bp, 0.0, 8)
        for lo_, hi_ in p.NUCLEOSOME_SITES_BP:
            ax.axvspan(lo_, hi_, color=PALETTE["light_grey"], alpha=0.5, lw=0, zorder=0)
        for x0, ls in ((387, ":"), (584, ":")):
            ax.axvline(x0, color=PALETTE["dark_grey"], lw=0.6, ls=ls, zorder=1)
        ax.axvline(1028 if SCALE == "norm1028" else 971, color="black", lw=0.6, ls="--")
        ax.bar(centers, sq / smax, width=5.0, color=PALETTE["grey"], lw=0, alpha=0.6)
        ax.plot(centers, base / smax, color=PALETTE["dark_grey"], lw=0.7, ls="--")
        ax.plot(centers, sig / smax, color=PALETTE["blue"], lw=0.9)
        for x in peaks:
            ax.plot(x, 1.02, "v", ms=4, color=PALETTE["red"], clip_on=False)
            ax.text(x, 1.06, f"{x:.0f}", color=PALETTE["red"], fontsize=6, ha="center", va="bottom")
            rows.append({"group": lab, "peak_bp": round(x, 0), "n_tethers": n_t})
        ax.text(0.01, 0.95, f"{lab}  ({n_t} tethers)", transform=ax.transAxes, fontsize=7, va="top")
        ax.set_ylim(0, 1.15)
        style_axes(ax, ylabel="sqrt (norm.)")
    style_axes(axes[-1], xlabel=("Normalised position (0-1028, per tether)" if SCALE == "norm1028"
                                 else f"Luz, model {MODEL} (bp from the 10 pN zero)"), ylabel="sqrt (norm.)")
    axes[-1].set_xlim(0, 1100)
    fig.suptitle(f"20261009, your {len(pooled)} selected tethers: "
                 f"{'0-1028 scale' if SCALE == 'norm1028' else 'Luz model ' + MODEL} peak map "
                 "(601 sites shaded, dyads dotted, 971 bp trunk end dashed)", fontsize=8)
    fig.tight_layout()
    for ext_ in ("pdf", "png"):
        fig.savefig(out / f"fig_{'norm1028' if SCALE == 'norm1028' else 'luz' + MODEL}_peakmap_pooled.{ext_}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out / f"{'norm1028' if SCALE == 'norm1028' else 'luz' + MODEL}_pooled_peaks.csv",
                              index=False)
    print("pooled peaks:", pd.DataFrame(rows).groupby("group").peak_bp.apply(list).to_dict())


def main() -> int:
    p.apply_style()
    sel = pd.read_csv(SEL_CSV)
    levels = sorted({10.0, 18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9, 46.0})
    colours = p.level_colours(levels)
    merged = fitz.open()
    for pdf in sorted(BACKUP.glob("*.pdf")):
        merged.insert_pdf(fitz.open(pdf))
    merged.save(p.DEST / f"WHZ_20261009_selected_{len(sel)}_FEC_pages.pdf", garbage=4, deflate=True)
    rows, toc, n, pooled = [], [], 0, []
    out_pdf = p.DEST / f"WHZ_20261009_selected_{len(sel)}_norm1028{SUFFIX}.pdf"
    with PdfPages(out_pdf) as pdf:
        for d, g in sel.groupby("data_id", sort=True):
            run = p.load(int(d))
            steps = p.segment_steps(run)
            trk = file_tracking(p.RAW / f"data_{int(d):03d}.h5")
            toc.append([1, f"data_{int(d):03d} ({len(g)})", n + 1])
            for r in g.sort_values("page_in_trimmed").itertuples():
                roi = r.bead[1:]
                tq = trk.get(roi, {})
                ext, rec = bead_scale(run, steps, roi, tq)
                title = f"20261009 data_{int(d):03d} {r.bead}   tracking {tq.get('status', '?')}"
                if np.isfinite(rec["end_nm"]):
                    rec["end_luz"] = float(model_luz(np.array([rec["end_nm"] - rec["zero_nm"]]),
                                                     np.array([rec["ref_force_pn"]]))[0])
                fig, peaks = draw(run, steps, roi, ext, rec, colours, title)
                if "luz_samples" in rec:
                    pooled.append((f"D{int(d):03d}_{r.bead}",) + rec["luz_samples"])
                pdf.savefig(fig, dpi=200, bbox_inches="tight")
                plt.close(fig)
                n += 1
                toc.append([2, f"{r.bead}  span {rec['span_nm']:.0f} nm", n])
                rows.append({"data_id": int(d), "bead": r.bead, "tracking": tq.get("status", "?"),
                             "t_loss_s": round(rec["t_loss"] - run["time_s"][0], 1) if np.isfinite(rec["t_loss"]) else np.nan,
                             "zero_10pN_max_nm": round(rec["zero_nm"], 1),
                             "ref_step": rec["ref_step"], "ref_force_pn": rec["ref_force_pn"],
                             "end_max_nm": round(rec["end_nm"], 1), "span_nm": round(rec["span_nm"], 1),
                             "skipped_steps": rec["skipped"],
                             "max_held_s": round(rec.get("max_held_s", np.nan), 2),
                             f"end_luz{MODEL}_bp": round(rec.get("end_luz", np.nan), 1),
                             ("peaks_norm1028" if SCALE == "norm1028" else f"luz{MODEL}_peaks_bp"):
                                 "; ".join(f"{x:.0f}" for x in peaks),
                             "drift": rec["drift"], "flags": " | ".join(rec["flags"])})
            print(f"data_{int(d):03d}: {len(g)}", flush=True)
    doc = fitz.open(out_pdf)
    doc.set_toc(toc)
    doc.saveIncr()
    tab = pd.DataFrame(rows)
    tab.to_csv(p.DEST / f"norm1028_selected_{len(sel)}{SUFFIX}.csv", index=False)
    pooled_peakmap(pooled)
    print(tab.ref_force_pn.value_counts().sort_index().to_dict())
    print("1028 from an earlier step:", int((tab.skipped_steps.fillna("") != "").sum()))
    print(tab.span_nm.describe().round(0).to_dict())
    print(f"-> {out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
