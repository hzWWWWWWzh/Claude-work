#!/usr/bin/env python3
"""20261004: the same tether at every maximum force - where does the fork dwell?

In this titration each field of view is pulled with a stepped-up Fmax
(18.0 -> 46.0 pN on the batch calibration the user gave on 2026-10-07:
1.7 / 1.0 / 0.7 mm = 18 / 30 / 46 pN). Every pull holds Fmax for ~10 s
(35-61 frames/s depending on the file; FOV02's holds are only 1-3 s).
That plateau is constant force, so the Luz samples on it are a dwell-time
histogram: where the fork sits, and for how long, at that force. The SOP branch
split (dShift/dt) leaves the plateau out - it is neither up nor down - so the
tables of process_20261004.py do not contain it; this script reads it back
from the raw H5 files.

Every plateau sample is normalised with the bead's single length scale and
the zero of its own file and cycle (the anchors in per_cycle_anchors.csv), so
within one tether the files are on the same bp axis and the ~11% uncertainty
of a borrowed scale is a common factor rather than a file-to-file error.

Stall position at a given Fmax = median Luz over the last 1 s of the plateau
(where the fork is after ~10 s at that force). A cycle counts ("used") only if
  - the hold lasted >= 9 s (FOV02 is excluded: its 1-3 s holds are not
    comparable),
  - the bead is still tracked at the end of the hold (a bead lost mid-plateau
    would otherwise report its position seconds earlier),
  - the fork is stationary in that last second (IQR <= 100 bp, the stall-window
    criterion of stall_positions_20261004.py), and
  - the value is within -150..1470 bp (a closed fork reads within ~50 bp of 0).
Tethers whose used stalls never exceed 50 bp at any force never unzipped; they
are reported, not mixed into the paired statistics.
Dwell times weight every sample by its own frame interval.

Statistics with incomplete ladders (no tether needs every force):
  - additive model stall = tether + force level (least squares), so each
    force's effect is estimated within tethers; 95% CI by resampling tethers;
  - each step between neighbouring forces of a field's ladder, on every tether
    that has both (Wilcoxon, Bonferroni over the steps tested);
  - per tether: its highest vs its lowest force (Wilcoxon), and Spearman's rho
    of stall against force for tethers with >= 3 forces.

Caveat built into the design: in every field Fmax rises with file number, so
a higher force is also a later pull of the same molecule. A fork that goes
further at higher Fmax may reflect force or the history of earlier pulls.
Within a file, cycle 2 (a repeat at the same force after a ~2 s rest) does
not go further than cycle 1.

Outputs (WHZ_20261004/tether_titration/):
  WHZ_20261004_tether_stitched_by_fmax.pdf   one page per tether: every file of
      the tether on one time axis in acquisition order (Luz and force, raw
      samples), each Fmax in its own colour, file boundaries and the analysed
      stall marked, and the plateau dwell histogram of each force side by side
      in the same colours (user's request, 2026-10-07)
  WHZ_20261004_tether_dwell_by_fmax.pdf   one page per tether, one row per file
      of the field: C1 and C2 plateau traces apart, and the dwell histogram
  fig_stall_by_fmax_all_tethers          every tether with >= 2 forces (lines join
      one tether) and the tether-adjusted mean per force with 95% CI
  fig_stall_histograms_by_fmax_all       stall histograms per force, all tethers
      that have that force
  trend_stats_all_tethers.csv, trend_sensitivity.csv
  fig_stall_tether_by_fmax_heatmap        tether x Fmax
  fig_stitched_FOVxx_Rxxx                three example tethers
  plateau_samples.csv.gz, stitched_luz_samples.csv.gz, plateau_stall_by_cycle.csv,
  plateau_stall_by_tether_level.csv, paired_stats_*.csv, heatmap_rows.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullLocator
from scipy.stats import wilcoxon
from whz_plot import PALETTE, apply_style, export_panel, raw_trace, style_axes

from process_tierA_sop import drift_corrected_nm
from process_20261004 import (
    CALIBRATION, CONSTRUCT_BP, DEST, FOVS, NUCLEOSOME_SITES_BP, RAW, normalise,
)
from unzip_raw_points import load_raw_run

OUT = DEST / "tether_titration"
PLATEAU_TOL_PN = 0.5        # plateau = force within this of the file's Fmax
MIN_PLATEAU_S = 2.0         # shorter holds are not plateaus at all (duration, not samples)
REQUIRED_PLATEAU_S = 9.0    # a stall is comparable only after the full ~10 s hold
END_WINDOW_S = 1.0          # stall = median Luz over the last second of the plateau
END_REACH_S = 0.1           # the kept samples must reach the plateau end this closely
MAX_END_IQR_BP = 100.0      # stationary fork (same rule as the 30-38 pN stall window)
STALL_GUARD_BP = (-150.0, 1.5 * CONSTRUCT_BP)
NEVER_UNZIPS_BP = 50.0
MOVING_SLOPE_BP_S = 20.0    # reported only: fork still advancing over the last 3 s
HIST_BIN_BP = 10.0
# A tether enters the comparison with any two forces (user, 2026-10-08: "not
# every force has to be there"); ladders are incomplete, so the statistics use
# all available pairs instead of complete cases.
MIN_LEVELS = 2
N_BOOT = 2000
LUZ_LIM = (-300.0, 1.15 * CONSTRUCT_BP)
CYC = {1: PALETTE["red"], 2: PALETTE["blue"]}
FILE_GAP_S = 8.0            # blank space between concatenated files
SITE1_ENTRY = NUCLEOSOME_SITES_BP[0][0]


def level_colours(levels: list[float]) -> dict[float, tuple]:
    """Fmax is ordered: one viridis step per level (pale end dropped), so the
    close middle levels stay distinct and one force has one colour everywhere."""
    cols = matplotlib.colormaps["viridis"](np.linspace(0, 0.85, len(levels)))
    return {f: tuple(c) for f, c in zip(levels, cols)}


def plateaus(run) -> dict[int, np.ndarray]:
    """Sample indices of each Fmax plateau, numbered as cycles in time order."""
    on = run.force_pn >= run.force_pn.max() - PLATEAU_TOL_PN
    idx = np.flatnonzero(on)
    if not idx.size:
        return {}
    segs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    t = run.time_s
    segs = [s for s in segs if t[s[-1]] - t[s[0]] >= MIN_PLATEAU_S]
    return {k + 1: s for k, s in enumerate(segs)}


def frame_dt(t: np.ndarray) -> np.ndarray:
    """Each sample's own frame interval (files run at 35-61 frames/s, unevenly)."""
    return np.gradient(t) if t.size > 1 else np.full(t.size, np.nan)


def usable_anchors() -> pd.DataFrame:
    a = pd.read_csv(DEST / "per_cycle_anchors.csv", low_memory=False)
    a["ref_cycles_disagree"] = a.ref_cycles_disagree.astype(str).eq("True")
    return a[a.anchor_valid & a.ref_kind.isin(["full", "partial"]) & ~a.ref_cycles_disagree]


def cycle_of_samples(run) -> np.ndarray:
    """Cycle of every sample: from the start of its up branch to the next one."""
    starts = sorted((float(run.time_s[seg[0]]), int(c)) for c, br, seg in run.branches if br == "up")
    out = np.full(run.time_s.size, starts[0][1] if starts else 1)
    for t0, c in starts:
        out[run.time_s >= t0] = c
    return out


def collect(anchors: pd.DataFrame):
    samples, stalls, plateau_log, whole, file_info = [], [], [], [], {}
    for fov, files in FOVS.items():
        for d in files:
            g = anchors[anchors.data_id == d]
            run = load_raw_run(RAW, d, CALIBRATION)
            pl = plateaus(run)
            fmax = float(run.force_pn.max())
            t0 = run.time_s[0]
            file_info[d] = {"fov": fov, "fmax_pn": round(fmax, 1), "time_s": run.time_s - t0,
                            "force_pn": run.force_pn, "normalised": not g.empty}
            plateau_log.append({"data_id": d, "fov": fov, "fmax_pn": round(fmax, 1),
                                "n_plateaus": len(pl),
                                "plateau_s": [round(float(run.time_s[s[-1]] - run.time_s[s[0]]), 2)
                                              for s in pl.values()],
                                "frames_per_s": round((len(run.time_s) - 1)
                                                      / (run.time_s[-1] - t0), 1)})
            if g.empty:
                continue
            cyc_of = cycle_of_samples(run)
            for r in g.itertuples():
                roi = r.roi_key.split("_R")[1]
                ext_all = drift_corrected_nm(run.z_um[roi], run.time_s)
                lost = r.tracking == "lost" and np.isfinite(r.t_loss_s)
                # the whole cycle - ramps and plateau - for the stitched trace
                m = cyc_of == int(r.cycle)
                _, luz_all = normalise(ext_all[m], run.force_pn[m], r.z_zero_raw_nm,
                                       r.raw_anchor_span_nm, r.zero_force_pn, r.full_force_pn)
                ok = np.isfinite(luz_all)
                if lost:
                    ok &= run.time_s[m] < r.t_loss_s
                if ok.any():
                    whole.append(pd.DataFrame({
                        "fov": fov, "bead": r.bead, "data_id": d, "fmax_pn": round(fmax, 1),
                        "cycle": int(r.cycle),
                        "t_file_s": np.round((run.time_s[m] - t0)[ok], 3),
                        "force_pn": np.round(run.force_pn[m][ok], 2),
                        "luz_bp": np.round(luz_all[ok], 1)}))
                seg = pl.get(int(r.cycle))
                if seg is None:
                    continue
                t = run.time_s[seg]
                _, luz = normalise(ext_all[seg], run.force_pn[seg], r.z_zero_raw_nm,
                                   r.raw_anchor_span_nm, r.zero_force_pn, r.full_force_pn)
                dt = frame_dt(t)
                keep = np.isfinite(luz)
                if lost:
                    keep &= t < r.t_loss_s
                if keep.sum() < 10:
                    continue
                hold_s = float(t[-1] - t[0])
                t_rel = t[keep] - t[0]
                reached_end = t_rel[-1] >= hold_s - END_REACH_S
                luz, dt = luz[keep], dt[keep]
                end = luz[t_rel >= hold_s - END_WINDOW_S] if reached_end else luz[:0]
                start = luz[t_rel <= 0.5]
                if end.size >= 3:
                    q1, q3 = np.percentile(end, [25, 75])
                    stall = float(np.median(end))
                else:
                    q1 = q3 = stall = np.nan
                last3 = t_rel >= hold_s - 3.0
                slope = (float(np.polyfit(t_rel[last3], luz[last3], 1)[0])
                         if reached_end and last3.sum() >= 10 else np.nan)
                why = ("hold < 9 s" if hold_s < REQUIRED_PLATEAU_S else
                       "lost before the end of the hold" if not reached_end else
                       "fork not stationary (IQR > 100 bp)" if not q3 - q1 <= MAX_END_IQR_BP else
                       "outside -150..1470 bp" if not (STALL_GUARD_BP[0] <= stall
                                                       <= STALL_GUARD_BP[1]) else "")
                stalls.append({"fov": fov, "bead": r.bead, "data_id": d, "fmax_pn": round(fmax, 1),
                               "cycle": int(r.cycle), "kind": r.ref_kind,
                               "ref_data_id": int(r.ref_data_id), "hold_s": round(hold_s, 2),
                               "tracked_s": round(float(t_rel[-1]), 2), "n_samples": int(keep.sum()),
                               "t_end_file_s": round(float(t[-1] - t0), 3),
                               "luz_start_bp": round(float(np.median(start)), 1) if start.size else np.nan,
                               "luz_end_bp": round(stall, 1), "end_iqr_bp": round(float(q3 - q1), 1),
                               "end_slope_bp_s": round(slope, 1),
                               "still_moving": bool(abs(slope) > MOVING_SLOPE_BP_S)
                               if np.isfinite(slope) else False,
                               "used": why == "", "why_not_used": why})
                samples.append(pd.DataFrame({"fov": fov, "bead": r.bead, "data_id": d,
                                             "fmax_pn": round(fmax, 1), "cycle": int(r.cycle),
                                             "t_plateau_s": np.round(t_rel, 3),
                                             "dt_s": np.round(dt, 5),
                                             "luz_bp": np.round(luz, 1)}))
            print(f"{fov} data_{d:03d} Fmax {fmax:4.1f} pN: {len(pl)} plateaus, "
                  f"{len(g)} anchored cycles", flush=True)
    pd.DataFrame(plateau_log).to_csv(OUT / "plateaus_by_file.csv", index=False)
    return (pd.concat(samples, ignore_index=True), pd.DataFrame(stalls),
            pd.concat(whole, ignore_index=True), file_info)


def tether_levels(st: pd.DataFrame) -> pd.DataFrame:
    """One stall per tether and force: median of its used cycles (both kept)."""
    u = st[st.used]
    tl = (u.groupby(["fov", "bead", "fmax_pn"], as_index=False)
           .agg(kind=("kind", "first"), luz_end_bp=("luz_end_bp", "median"),
                n_cycles=("cycle", "nunique"), data_id=("data_id", "first")))
    by_c = (u.pivot_table(index=["fov", "bead", "fmax_pn"], columns="cycle", values="luz_end_bp")
             .rename(columns={1: "c1_bp", 2: "c2_bp"}).reset_index())
    by_c.columns.name = None
    tl = tl.merge(by_c, how="left")
    reach = tl.groupby(["fov", "bead"]).luz_end_bp.transform("max")
    tl["never_unzips"] = reach <= NEVER_UNZIPS_BP
    return tl


def shade_sites(ax, horizontal: bool):
    for lo, hi in NUCLEOSOME_SITES_BP:
        (ax.axhspan if horizontal else ax.axvspan)(lo, hi, color=PALETTE["light_grey"],
                                                   alpha=0.5, lw=0, zorder=0)


def dwell_hist(luz: pd.Series, dt: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(LUZ_LIM[0], LUZ_LIM[1] + HIST_BIN_BP, HIST_BIN_BP)
    h, _ = np.histogram(luz.clip(*LUZ_LIM), bins=edges, weights=dt)
    return h, edges


def tether_pages(samples, st, file_info, keep, pdf_path: Path) -> int:
    """Per tether, one row per file of its field: C1 trace, C2 trace, dwell histogram."""
    n = 0
    with PdfPages(pdf_path) as pdf:
        for fov, bead in keep:
            s = samples[(samples.fov == fov) & (samples.bead == bead)]
            c = st[(st.fov == fov) & (st.bead == bead)]
            files = [d for d in FOVS[fov] if d in file_info]
            ref = int(c.ref_data_id.iloc[0])
            fig, axes = plt.subplots(len(files), 3, figsize=(7.09, 1.35 * len(files) + 0.6),
                                     gridspec_kw={"width_ratios": [1.5, 1.5, 1]},
                                     sharey=True, squeeze=False)
            hmax = 0.0
            for k, d in enumerate(files):
                f = file_info[d]["fmax_pn"]
                for j, ax in enumerate(axes[k]):
                    shade_sites(ax, True)
                    ax.set_ylim(*LUZ_LIM)
                for cyc in (1, 2):
                    ax = axes[k, cyc - 1]
                    g = s[(s.data_id == d) & (s.cycle == cyc)]
                    cc = c[(c.data_id == d) & (c.cycle == cyc)]
                    if g.empty:
                        ax.text(0.5, 0.5, "not normalised" if not len(cc) else "no plateau",
                                transform=ax.transAxes, ha="center", va="center",
                                color=PALETTE["dark_grey"])
                        continue
                    hold = float(cc.hold_s.iloc[0]) if len(cc) else g.t_plateau_s.max()
                    ax.axvspan(hold - END_WINDOW_S, hold, color=PALETTE["light_grey"],
                               alpha=0.35, lw=0, zorder=0)
                    raw_trace(ax, g.t_plateau_s, g.luz_bp, color=CYC[cyc], marker_area=1.5,
                              alpha=0.6, edge_width=0.3, rasterized=True)
                    h, edges = dwell_hist(g.luz_bp, g.dt_s)
                    axes[k, 2].stairs(h, edges, orientation="horizontal", color=CYC[cyc], lw=0.9)
                    hmax = max(hmax, h.max())
                    if len(cc) and np.isfinite(cc.luz_end_bp.iloc[0]):
                        axes[k, 2].axhline(cc.luz_end_bp.iloc[0], color=CYC[cyc], lw=0.7,
                                           ls="--" if cc.used.iloc[0] else ":")
                axes[k, 0].set_title(f"data_{d:03d}  {f:.1f} pN"
                                     + ("  [length ref]" if d == ref else ""), loc="left")
                last = k == len(files) - 1
                style_axes(axes[k, 0], xlabel="C1 time on plateau (s)" if last else None,
                           ylabel="Luz (bp)")
                style_axes(axes[k, 1], xlabel="C2 time on plateau (s)" if last else None)
                style_axes(axes[k, 2], xlabel="Dwell per 10 bp (s)" if last else None)
            for ax in axes[:, 2]:
                ax.set_xlim(0, hmax * 1.05 if hmax else 1)
            for ax in axes[:, 0:2].flat:
                ax.set_xlim(-0.3, 10.5)
            fig.suptitle(f"20261004 {fov} {bead} ({c.kind.iloc[0]} unzip). Shaded: last 1 s; "
                         "dashed = stall used, dotted = not used", y=1.0)
            fig.tight_layout(h_pad=0.6, w_pad=0.5)
            pdf.savefig(fig, dpi=200, bbox_inches="tight")
            plt.close(fig)
            n += 1
    return n


def draw_stitched(fig, bead_whole, bead_plateau, bead_stalls, bead_levels, fov, file_info, cols):
    """Every file of one tether on one time axis; per-force dwell side by side."""
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], width_ratios=[3.6, 1.4],
                          hspace=0.08, wspace=0.05)
    ax_l = fig.add_subplot(gs[0, 0])
    # not sharex: the Luz panel's ticks label the files (top), the force panel's the time
    ax_f = fig.add_subplot(gs[1, 0])
    ax_h = fig.add_subplot(gs[0, 1], sharey=ax_l)
    shade_sites(ax_l, True)
    shade_sites(ax_h, True)
    offset = 0.0
    centres, labels = [], []
    files = [d for d in FOVS[fov] if d in file_info]
    for d in files:
        info = file_info[d]
        f = info["fmax_pn"]
        raw_trace(ax_f, info["time_s"] + offset, info["force_pn"], color=cols[f], marker_area=0.4,
                  alpha=0.8, edge_width=0.1, filled=True, rasterized=True)
        g = bead_whole[bead_whole.data_id == d]
        if len(g):
            raw_trace(ax_l, g.t_file_s + offset, g.luz_bp, color=cols[f], marker_area=0.8,
                      alpha=0.55, edge_width=0.15, filled=True, rasterized=True)
        for r in bead_stalls[bead_stalls.data_id == d].itertuples():
            if np.isfinite(r.luz_end_bp):
                ax_l.hlines(r.luz_end_bp, offset + r.t_end_file_s - END_WINDOW_S,
                            offset + r.t_end_file_s, color="black" if r.used else PALETTE["grey"],
                            lw=1.6, zorder=4)
        centres.append(offset + info["time_s"][-1] / 2)
        labels.append(f"{f:.1f}" + ("" if len(g) else "*"))
        end = offset + info["time_s"][-1]
        offset = end + FILE_GAP_S
        if d != files[-1]:
            for ax in (ax_l, ax_f):
                ax.axvline(end + FILE_GAP_S / 2, color=PALETTE["light_grey"], lw=0.6, zorder=0)
    # dwell at each force, side by side on a categorical axis that shares Luz
    levels = sorted(bead_plateau.fmax_pn.unique())
    hists = {f: dwell_hist(g.luz_bp, g.dt_s) for f, g in bead_plateau.groupby("fmax_pn")}
    hmax = max((h.max() for h, _ in hists.values()), default=1.0)
    for k, f in enumerate(levels):
        h, edges = hists[f]
        mids = (edges[:-1] + edges[1:]) / 2
        ax_h.fill_betweenx(mids, k, k + 0.85 * h / hmax, step="mid", color=cols[f],
                           alpha=0.8, lw=0)
        # each cycle's own stall: the median of two different cycles can be a
        # position the fork never occupied
        for r in bead_stalls[(bead_stalls.fmax_pn == f)].itertuples():
            if np.isfinite(r.luz_end_bp):
                ax_h.hlines(r.luz_end_bp, k, k + 0.85, lw=1.2, zorder=4,
                            color="black" if r.used else PALETTE["grey"])
    ax_h.set_xlim(-0.15, len(levels))
    ax_h.set_xticks(np.arange(len(levels)) + 0.42, [f"{f:.1f}" for f in levels], rotation=90)
    ax_h.xaxis.set_minor_locator(NullLocator())
    ax_l.set_ylim(*LUZ_LIM)
    for ax in (ax_l, ax_f):
        ax.set_xlim(-1, offset - FILE_GAP_S + 1)
    ax_f.set_ylim(0, max(cols) * 1.1)
    style_axes(ax_l, ylabel="Luz (bp)")
    style_axes(ax_f, xlabel="Time, files in acquisition order (s)", ylabel="Force (pN)")
    style_axes(ax_h, xlabel="Fmax (pN)")
    # file labels (Fmax, * = not normalised) as the Luz panel's own top tick labels
    ax_l.set_xticks(centres, labels)
    ax_l.xaxis.set_minor_locator(NullLocator())
    ax_l.tick_params(axis="x", labelbottom=False, labeltop=True)
    ax_h.tick_params(labelleft=False)
    return ax_l


def stitched_pages(whole, samples, st, tl, file_info, keep, cols, pdf_path, examples=()):
    def one(fig, fov, bead):
        sel = lambda df: df[(df.fov == fov) & (df.bead == bead)]
        draw_stitched(fig, sel(whole), sel(samples), sel(st), sel(tl), fov, file_info, cols)

    n = 0
    with PdfPages(pdf_path) as pdf:
        for fov, bead in keep:
            kind = st[(st.fov == fov) & (st.bead == bead)].kind.iloc[0]
            fig = plt.figure(figsize=(7.09, 3.6))
            one(fig, fov, bead)
            fig.suptitle(f"20261004 {fov} {bead} ({kind} unzip), files in acquisition order "
                         "labelled by Fmax (pN), * = not normalised; black bar = stall used, "
                         "grey = not used; right: dwell at each Fmax", y=1.04, fontsize=8)
            pdf.savefig(fig, dpi=200, bbox_inches="tight")
            plt.close(fig)
            n += 1
    for fov, bead in examples:
        fig = plt.figure(figsize=(7.09, 3.4))
        one(fig, fov, bead)
        # stacked panels share axes, so the per-axis label check is the only one off
        export_panel(fig, OUT / f"fig_stitched_{fov}_{bead}", dpi=600, require_axis_labels=False)
    return n


def analysis_set(tl: pd.DataFrame) -> pd.DataFrame:
    """Unzipping tethers with a used stall at >= MIN_LEVELS forces (long format)."""
    t = tl[~tl.never_unzips].copy()
    t["tether"] = t.fov + "_" + t.bead
    n = t.groupby("tether").fmax_pn.transform("nunique")
    return t[n >= MIN_LEVELS].reset_index(drop=True)


def additive_fit(t: pd.DataFrame, levels: list[float]) -> dict[float, float]:
    """stall = tether + level, least squares; returns the tether-averaged
    prediction at every level (the level effect is estimated within tethers)."""
    tethers = sorted(t.tether.unique())
    ti = {k: i for i, k in enumerate(tethers)}
    present = [f for f in levels if f in set(t.fmax_pn)]
    li = {f: j for j, f in enumerate(present)}
    x = np.zeros((len(t), len(tethers) + len(present) - 1))
    for r, (k, f) in enumerate(zip(t.tether, t.fmax_pn)):
        x[r, ti[k]] = 1.0
        if li[f]:
            x[r, len(tethers) + li[f] - 1] = 1.0
    beta = np.linalg.lstsq(x, t.luz_end_bp.to_numpy(float), rcond=None)[0]
    base = beta[:len(tethers)].mean()
    eff = np.r_[0.0, beta[len(tethers):]]
    return {f: float(base + eff[li[f]]) if f in li else np.nan for f in levels}


def bootstrap_fit(t: pd.DataFrame, levels: list[float]) -> pd.DataFrame:
    """95% CI of the adjusted means and of their shift from the lowest force,
    resampling whole tethers."""
    rng = np.random.default_rng(0)
    groups = {k: g for k, g in t.groupby("tether")}
    keys = list(groups)
    boots = []
    for b in range(N_BOOT):
        pick = rng.choice(len(keys), len(keys), replace=True)
        bt = pd.concat([groups[keys[i]].assign(tether=f"{keys[i]}#{j}")
                        for j, i in enumerate(pick)], ignore_index=True)
        boots.append(additive_fit(bt, levels))
    bt = pd.DataFrame(boots)
    shift = bt.sub(bt[levels[0]], axis=0)
    fit = additive_fit(t, levels)
    rows = []
    for f in levels:
        rows.append({"fmax_pn": f, "adjusted_mean_bp": round(fit[f], 0),
                     "ci_lo_bp": round(bt[f].quantile(0.025), 0),
                     "ci_hi_bp": round(bt[f].quantile(0.975), 0),
                     "shift_from_lowest_bp": round(fit[f] - fit[levels[0]], 0),
                     "shift_ci_lo_bp": round(shift[f].quantile(0.025), 0),
                     "shift_ci_hi_bp": round(shift[f].quantile(0.975), 0)})
    return pd.DataFrame(rows)


def step_tests(t: pd.DataFrame, ladders: dict[str, list[float]]) -> pd.DataFrame:
    """Neighbouring forces of each field's ladder, on every tether that has both."""
    wide = t.pivot_table(index=["fov", "tether"], columns="fmax_pn", values="luz_end_bp")
    diffs = {}
    for (fov, _), row in wide.iterrows():
        lad = ladders[fov]
        for a, b in zip(lad[:-1], lad[1:]):
            if a in row and b in row and np.isfinite(row[a]) and np.isfinite(row[b]):
                diffs.setdefault((a, b), []).append(row[b] - row[a])
    rows = []
    for (a, b), d in sorted(diffs.items()):
        d = np.array(d)
        p = wilcoxon(d).pvalue if len(d) >= 5 and (d != 0).any() else np.nan
        rows.append({"comparison": f"{a}->{b}", "n": len(d), "median_bp": round(np.median(d), 0),
                     "q25_bp": round(np.percentile(d, 25), 0),
                     "q75_bp": round(np.percentile(d, 75), 0),
                     "frac_up": round((d > 0).mean(), 2), "wilcoxon_p": p})
    out = pd.DataFrame(rows)
    tested = out.wilcoxon_p.notna().sum()
    out["p_bonferroni"] = (out.wilcoxon_p * tested).clip(upper=1.0)
    return out


def trend_tests(t: pd.DataFrame) -> pd.DataFrame:
    from scipy.stats import spearmanr
    lo_hi, rhos = [], []
    for _, g in t.groupby("tether"):
        g = g.sort_values("fmax_pn")
        lo_hi.append(g.luz_end_bp.iloc[-1] - g.luz_end_bp.iloc[0])
        if len(g) >= 3:
            rhos.append(spearmanr(g.fmax_pn, g.luz_end_bp)[0])
    lo_hi, rhos = np.array(lo_hi), np.array(rhos)
    rhos = rhos[np.isfinite(rhos)]
    return pd.DataFrame([
        {"comparison": "highest - lowest force of each tether", "n": len(lo_hi),
         "median_bp": round(np.median(lo_hi), 0), "q25_bp": round(np.percentile(lo_hi, 25), 0),
         "q75_bp": round(np.percentile(lo_hi, 75), 0), "frac_up": round((lo_hi > 0).mean(), 2),
         "wilcoxon_p": wilcoxon(lo_hi).pvalue},
        {"comparison": "Spearman rho, stall vs force (tethers with >= 3 forces)", "n": len(rhos),
         "median_bp": round(np.median(rhos), 2), "q25_bp": round(np.percentile(rhos, 25), 2),
         "q75_bp": round(np.percentile(rhos, 75), 2), "frac_up": round((rhos > 0).mean(), 2),
         "wilcoxon_p": wilcoxon(rhos).pvalue}])


def level_table(t: pd.DataFrame, levels: list[float], boot: pd.DataFrame) -> pd.DataFrame:
    inside = lambda v: any(lo <= v <= hi for lo, hi in NUCLEOSOME_SITES_BP)
    rows = []
    for f in levels:
        x = t[t.fmax_pn == f].luz_end_bp
        rows.append({"fmax_pn": f, "n_tethers": len(x), "median_bp": round(x.median(), 0),
                     "q25_bp": round(x.quantile(.25), 0), "q75_bp": round(x.quantile(.75), 0),
                     "frac_past_314": round((x > SITE1_ENTRY).mean(), 2),
                     "frac_inside_601": round(x.apply(inside).mean(), 2)})
    return pd.DataFrame(rows).merge(boot, on="fmax_pn")


def all_tether_figures(t: pd.DataFrame, levels: list[float], boot: pd.DataFrame):
    """Every tether with >= 2 forces, lines joining its own forces; black =
    tether-adjusted mean and 95% CI. Plus stall histograms per force."""
    pos = {f: k for k, f in enumerate(levels)}
    fig, ax = plt.subplots(figsize=(3.39, 2.8))
    shade_sites(ax, True)
    # one small sideways offset per tether, the same at every force, so a line
    # and its points stay together
    tethers = sorted(t.tether.unique())
    jit = dict(zip(tethers, np.random.default_rng(0).uniform(-0.15, 0.15, len(tethers))))
    for k, g in t.groupby("tether"):
        g = g.sort_values("fmax_pn")
        xs = g.fmax_pn.map(pos) + jit[k]
        ax.plot(xs, g.luz_end_bp, color=PALETTE["light_grey"], lw=0.5, zorder=1)
    raw_trace(ax, t.fmax_pn.map(pos) + t.tether.map(jit), t.luz_end_bp,
              color=PALETTE["dark_grey"], filled=True, marker_area=6, alpha=0.75)
    b = boot.set_index("fmax_pn").loc[levels]
    ax.errorbar(range(len(levels)), b.adjusted_mean_bp,
                yerr=[b.adjusted_mean_bp - b.ci_lo_bp, b.ci_hi_bp - b.adjusted_mean_bp],
                fmt="o", color="black", ms=4, capsize=2, zorder=5)
    ax.set_xticks(range(len(levels)), [f"{f:.1f}" for f in levels], rotation=90)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(-0.5, len(levels) - 0.5)
    ax.set_ylim(*LUZ_LIM)
    style_axes(ax, xlabel="Fmax (pN)", ylabel="Stall at end of plateau (bp)")
    export_panel(fig, OUT / "fig_stall_by_fmax_all_tethers", dpi=600)

    bins = np.arange(STALL_GUARD_BP[0], CONSTRUCT_BP + 151, 30)
    fig, axes = plt.subplots(len(levels), 1, figsize=(3.39, 0.85 * len(levels) + 0.6),
                             sharex=True, squeeze=False)
    ymax = 0
    for ax, f in zip(axes[:, 0], levels):
        x = t[t.fmax_pn == f].luz_end_bp
        shade_sites(ax, False)
        h, _, _ = ax.hist(x, bins=bins, histtype="step", color="black", lw=1.0)
        ymax = max(ymax, h.max())
        ax.axvline(x.median(), color=PALETTE["dark_grey"], lw=0.9, ls="--")
        ax.text(1.02, 0.5, f"{f:.1f}\nn={len(x)}", transform=ax.transAxes, ha="left",
                va="center", fontsize=8)
    for k, ax in enumerate(axes[:, 0]):
        ax.set_ylim(0, ymax * 1.15)
        ax.yaxis.set_minor_locator(NullLocator())
        style_axes(ax, xlabel="Stall at end of plateau (bp)" if k == len(levels) - 1 else None)
        ax.set_xlim(bins[0], bins[-1])
    fig.supylabel("Tethers (count)")
    fig.text(0.985, 0.5, "Fmax (pN)", rotation=270, va="center", ha="right")
    fig.tight_layout(h_pad=0.15, rect=(0, 0, 0.9, 1))
    # stacked panels share one x label and one y label; the per-axis label
    # check is the only validation switched off
    export_panel(fig, OUT / "fig_stall_histograms_by_fmax_all", dpi=600,
                 require_axis_labels=False)


def heatmap(tl: pd.DataFrame, levels: list[float]):
    t = tl[~tl.never_unzips]
    wide = t.pivot_table(index=["fov", "bead"], columns="fmax_pn", values="luz_end_bp")
    wide = wide[[f for f in levels if f in wide.columns]]
    wide = wide[wide.notna().sum(axis=1) >= MIN_LEVELS]
    wide = wide.loc[wide.max(axis=1).sort_values().index]
    wide.reset_index().to_csv(OUT / "heatmap_rows.csv", index=False)
    fig, ax = plt.subplots(figsize=(3.39, 4.6))
    # cividis: viridis already means Fmax in the stitched pages
    cmap = matplotlib.colormaps["cividis"].copy()
    cmap.set_bad("white")
    im = ax.imshow(np.ma.masked_invalid(wide.values), aspect="auto", cmap=cmap,
                   vmin=0, vmax=CONSTRUCT_BP, interpolation="nearest")
    ax.set_xticks(range(len(wide.columns)), [f"{c:.1f}" for c in wide.columns], rotation=90)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticks([])
    style_axes(ax, xlabel="Fmax (pN)", ylabel=f"Tethers (n = {len(wide)})")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, extend="both")
    cb.set_label("Stall at end of plateau (bp)")
    cb.ax.tick_params(direction="in")
    cb.outline.set_visible(False)
    for sp in cb.ax.spines.values():
        sp.set_visible(True)
    for lo, hi in NUCLEOSOME_SITES_BP:
        cb.ax.axhspan(lo, hi, facecolor="none", edgecolor="white", lw=0.6)
    # the colour bar carries its label on its long axis; the per-axis label
    # check is the only validation switched off
    export_panel(fig, OUT / "fig_stall_tether_by_fmax_heatmap", dpi=600,
                 require_axis_labels=False)
    return len(wide)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    apply_style()
    anchors = usable_anchors()
    samples, st, whole, file_info = collect(anchors)
    whole.to_csv(OUT / "stitched_luz_samples.csv.gz", index=False)
    samples.to_csv(OUT / "plateau_samples.csv.gz", index=False)
    st.to_csv(OUT / "plateau_stall_by_cycle.csv", index=False)
    print(f"\nplateau cycles: {len(st)}; used: {int(st.used.sum())}; not used: "
          f"{st[~st.used].why_not_used.value_counts().to_dict()}; used but still moving "
          f"(> {MOVING_SLOPE_BP_S:.0f} bp/s over the last 3 s): "
          f"{int((st.used & st.still_moving).sum())}")
    tl = tether_levels(st)
    tl.to_csv(OUT / "plateau_stall_by_tether_level.csv", index=False)
    never = tl[tl.never_unzips].drop_duplicates(["fov", "bead"])
    print(f"tethers that never unzip (no used stall > {NEVER_UNZIPS_BP:.0f} bp): "
          f"{len(never)} {list(zip(never.fov, never.bead))}")

    nlev = st[st.used].groupby(["fov", "bead"]).fmax_pn.nunique()
    keep = sorted(nlev[nlev >= MIN_LEVELS].index)
    # every Fmax drawn on a page gets its colour: all files of the fields shown
    levels = sorted({file_info[d]["fmax_pn"] for fov in {k[0] for k in keep}
                     for d in FOVS[fov] if d in file_info})
    cols = level_colours(levels)
    print(f"tether pages (>= {MIN_LEVELS} forces): "
          f"{tether_pages(samples, st, file_info, keep, OUT / 'WHZ_20261004_tether_dwell_by_fmax.pdf')}")
    # examples: among unzipping tethers with the most forces, the one whose stall
    # moves most, the one that moves least, and a full unzip
    span = (tl[~tl.never_unzips].groupby(["fov", "bead"])
            .agg(n=("fmax_pn", "nunique"), rng=("luz_end_bp", lambda x: x.max() - x.min()),
                 kind=("kind", "first")))
    top = span[span.n == span.n.max()]
    ex = [top.rng.idxmax(), top.rng.idxmin()]
    full = span[(span.kind == "full") & (span.n >= MIN_LEVELS)]
    if len(full):
        ex.append(full.n.idxmax())
    n_st = stitched_pages(whole, samples, st, tl, file_info, keep, cols,
                          OUT / "WHZ_20261004_tether_stitched_by_fmax.pdf", examples=ex)
    print(f"stitched pages: {n_st}; example panels: {ex}")

    # ---- every tether with >= 2 forces; ladders need not be complete ----
    t = analysis_set(tl)
    lv = sorted(t.fmax_pn.unique())
    ladders = {fov: sorted({file_info[d]["fmax_pn"] for d in FOVS[fov] if d in file_info})
               for fov in t.fov.unique()}
    boot = bootstrap_fit(t, lv)
    levels_tab = level_table(t, lv, boot)
    steps = step_tests(t, ladders)
    trend = trend_tests(t)
    all_tether_figures(t, lv, boot)
    pd.concat([levels_tab.assign(section="force level"),
               steps.assign(section="neighbouring forces, same tether"),
               trend.assign(section="per-tether trend")], ignore_index=True).to_csv(
        OUT / "trend_stats_all_tethers.csv", index=False)
    print(f"\nall tethers with >= {MIN_LEVELS} forces: {t.tether.nunique()} "
          f"{t.drop_duplicates('tether').fov.value_counts().to_dict()}")
    print(levels_tab.to_string(index=False))
    print(steps.to_string(index=False))
    print(trend.to_string(index=False))

    # sensitivity: the same adjusted shift lowest -> highest force under other choices
    sens = []
    def shift(tt, label):
        if tt.tether.nunique() < 5:
            return
        b = bootstrap_fit(tt, sorted(tt.fmax_pn.unique()))
        top = b.iloc[-1]
        sens.append({"variant": label, "n_tethers": tt.tether.nunique(),
                     "levels": f"{b.fmax_pn.iloc[0]}-{top.fmax_pn}",
                     "shift_bp": top.shift_from_lowest_bp,
                     "ci_lo_bp": top.shift_ci_lo_bp, "ci_hi_bp": top.shift_ci_hi_bp})
    shift(t, "default (>= 2 forces)")
    loose = st.used | st.why_not_used.eq("fork not stationary (IQR > 100 bp)")
    shift(analysis_set(tether_levels(st.assign(used=loose))), "IQR rule off")
    shift(t[t.groupby("tether").fmax_pn.transform("nunique") >= 3], ">= 3 forces")
    # the longest ladder that most fields share (FOV06 ends at 42.8 instead of 46.0)
    lad_counts = pd.Series([tuple(v) for v in ladders.values()]).value_counts()
    ful = max(lad_counts.index, key=lambda lad: (len(lad), lad_counts[lad]))
    comp = t.groupby("tether").fmax_pn.apply(lambda x: set(ful) <= set(x))
    shift(t[t.tether.isin(comp[comp].index)], "complete ladder only")
    shift(t[t.kind == "partial"], "partial beads only")
    sens = pd.DataFrame(sens)
    sens.to_csv(OUT / "trend_sensitivity.csv", index=False)
    print(sens.to_string(index=False))
    print(f"\nheatmap tethers: {heatmap(tl, levels)}")

    c12 = st[st.used].pivot_table(index=["fov", "bead", "data_id"], columns="cycle",
                                  values="luz_end_bp").dropna()
    d = c12[2] - c12[1]
    print(f"\nwithin a file, C2 - C1 at the same force: n={len(d)}, median {d.median():.0f} bp "
          f"(IQR {d.quantile(.25):.0f}..{d.quantile(.75):.0f}), C2 > C1 in {(d > 0).mean():.0%}, "
          f"Wilcoxon p={wilcoxon(d).pvalue:.3f}")
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
