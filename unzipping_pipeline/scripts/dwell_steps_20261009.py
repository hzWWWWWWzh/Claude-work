#!/usr/bin/env python3
"""20261009: long stalls and the brief pauses between them, inside every plateau.

User, 2026-10-11: between two main (long) stalls there are brief jumps that
may be the nucleosome's outer and inner turns; analyse them. And: once the
force exceeds the critical force for one nucleosome, why does the fork still
wait in front of the second?

Signal. Each 10 s plateau is at constant force, so Luz is a linear function of
extension there and steps can be found on the extension itself, at the full
frame rate (30-41 frames/s, no median/mean smoothing, which would smear a
0.1 s dwell). Extension = raw Z minus the rest-dwell drift line, zero = the
maximum extension at 10 pN in step 1 (both as in fec_all_20261009.py). Luz is
the dWLC model with no per-bead length scale.

Step finding: binary segmentation on the extension, each split accepted when
it lowers the squared error by more than PENALTY x sigma^2 x ln(n) (BIC-like;
sigma from the MAD of first differences / sqrt 2 of that plateau); segments
shorter than MIN_FRAMES are not allowed; neighbours closer than MERGE_BP are
merged. A dwell is LONG if it lasts >= LONG_S, BRIEF if MIN_FRAMES frames to
LONG_S. Detection limit: ~0.1 s (4 frames at 40 frames/s).

Tethers: a valid 10 pN zero, tracked (lost beads only before the loss). The
file-level LUT level jump is recorded, not used: all 750-880 bp "stalls" of the
first analysis came from beads with 0.5-1.5 um jumps, but so do clean tethers
(data_015 R006, 1.14 um). Instead every PLATEAU is gated locally: if any fitted
level lies outside PLATEAU_GUARD_BP the bead sat on another LUT level during
that hold, and the whole plateau is dropped. Pass --selection CSV
(data_id, bead) to analyse the user's own pick instead.

Time resolution follows the camera, which slows at high force: 100 frames/s up
to ~28 pN, 50 at 30-34.6 pN, 25 at 37.1-39.9 pN (data_015), so MIN_FRAMES = 4
means 0.04 s, 0.08 s and 0.16 s. Frame noise 20-40 nm, i.e. 27-53 bp.

Nucleosome map (both 601 copies forward, tandem): N1 314-460 (dyad 387),
N2 511-657 (dyad 584); outer / inner turn boundaries at dyad -+40 bp.

    ~/.venvs/trace_editor_bm/bin/python dwell_steps_20261009.py [--selection picks.csv]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_20261009 as p  # noqa: E402
import fec_all_20261009 as F  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402
from reprocess_26_core_sop import HANDLE_BP, delta_nm_per_bp, x_ds  # noqa: E402

OUT = p.DEST / "dwell_steps"
PENALTY = 3.0
MIN_FRAMES = 4
MERGE_BP = 12.0
LONG_S = 1.0
PLATEAU_GUARD_BP = (-150.0, 1150.0)
N1, N2 = (314.0, 460.0), (511.0, 657.0)
DYADS = (387.0, 584.0)
TURN_EDGES = (347.0, 427.0, 544.0, 624.0)       # dyad -+ 40 bp
LANDMARKS = {"N1 entry": 314, "N1 outer|inner": 347, "N1 dyad": 387, "N1 inner|outer": 427,
             "N1 exit": 460, "N2 entry": 511, "N2 outer|inner": 544, "N2 dyad": 584,
             "N2 inner|outer": 624, "N2 exit": 657, "open": 980}


def segment(y: np.ndarray, sigma: float) -> list[tuple[int, int]]:
    """Binary segmentation of a piecewise-constant signal; returns [start, end)."""
    n = y.size
    cs, cs2 = np.concatenate([[0], np.cumsum(y)]), np.concatenate([[0], np.cumsum(y * y)])

    def sse(a, b):
        m = b - a
        s = cs[b] - cs[a]
        return (cs2[b] - cs2[a]) - s * s / m

    pen = PENALTY * sigma * sigma * np.log(max(n, 2))
    todo, done = [(0, n)], []
    while todo:
        a, b = todo.pop()
        if b - a < 2 * MIN_FRAMES:
            done.append((a, b)); continue
        ks = np.arange(a + MIN_FRAMES, b - MIN_FRAMES + 1)
        m1, m2 = ks - a, b - ks
        s1, s2 = cs[ks] - cs[a], cs[b] - cs[ks]
        q1, q2 = cs2[ks] - cs2[a], cs2[b] - cs2[ks]
        split = (q1 - s1 * s1 / m1) + (q2 - s2 * s2 / m2)
        k = int(ks[np.argmin(split)])
        if sse(a, b) - split.min() > pen:
            todo += [(a, k), (k, b)]
        else:
            done.append((a, b))
    return sorted(done)


def merge(segs, y, min_gap):
    out = []
    for a, b in segs:
        if out and abs(np.mean(y[a:b]) - np.mean(y[out[-1][0]:out[-1][1]])) < min_gap:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def tethers(selection: Path | None) -> pd.DataFrame | None:
    return pd.read_csv(selection) if selection else None


def plateau_dwells(d: int, picks: pd.DataFrame | None):
    run = p.load(d)
    steps = p.segment_steps(run)
    trk = file_tracking(p.RAW / f"data_{d:03d}.h5")
    t_abs = run["time_s"]
    t = t_abs - t_abs[0]
    rows, traces = [], {}
    for roi in run["rois"]:
        bead = f"R{roi}"
        if picks is not None and not ((picks.data_id == d) & (picks.bead == bead)).any():
            continue
        q = trk.get(roi, {})
        if picks is None and q.get("status") not in ("tracked", "lost"):
            continue
        t_loss = (q["t_loss_s"] - t_abs[0]) if q.get("status") == "lost" else np.inf
        z_raw = run["z_um"][roi] * 1000.0
        z_s = F.smoothed_nm(run["z_um"][roi])
        run["z_um_raw"] = run["z_um"][roi]
        ext_s, note = F.rest_drift(run, steps, z_s)
        drift = z_s - ext_s                                 # the drift line itself
        ext_raw = z_raw - drift
        pl1 = steps[0]["plateau"]
        zero = F.max_extension(t[pl1], ext_s[pl1])
        if not np.isfinite(zero):
            continue
        for st in steps:
            pl = st["plateau"]
            keep = (t[pl] < t_loss) & np.isfinite(ext_raw[pl])
            if keep.sum() < 3 * MIN_FRAMES:
                continue
            idx = pl[keep]
            fF = st["force_pn"]
            dx = float(delta_nm_per_bp(np.array([fF]))[0])
            hs = float(HANDLE_BP * (x_ds(np.array([fF])) - x_ds(np.array([F.ZERO_F_PN])))[0])
            luz = (ext_raw[idx] - zero - hs) / dx
            sig = float(np.median(np.abs(np.diff(luz))) / 0.6745 / np.sqrt(2)) or 1.0
            segs = merge(segment(luz, sig), luz, MERGE_BP)
            levels = [float(np.mean(luz[a:b])) for a, b in segs]
            if min(levels) < PLATEAU_GUARD_BP[0] or max(levels) > PLATEAU_GUARD_BP[1]:
                continue                                   # another LUT level in this hold
            tt = t[idx]
            frame = float(np.median(np.diff(tt)))
            for k, (a, b) in enumerate(segs):
                dur = tt[b - 1] - tt[a] + frame
                rows.append({"data_id": d, "bead": bead, "step": st["step"], "force_pn": round(fF, 1),
                             "k": k, "n_segments": len(segs), "t_start_s": round(tt[a] - tt[0], 3),
                             "duration_s": round(dur, 3), "level_bp": round(float(np.mean(luz[a:b])), 1),
                             "level_sem_bp": round(float(np.std(luz[a:b]) / np.sqrt(b - a)), 1),
                             "sigma_bp": round(sig, 1), "frames": b - a,
                             "censored": a == 0 or b == len(luz),
                             "fps": round(1.0 / frame, 0), "level_jump_um": q.get("level_jump_um", np.nan)})
            traces[(d, bead, st["step"])] = (tt - tt[0], luz, segs)
    return pd.DataFrame(rows), traces


def classify(dw: pd.DataFrame) -> pd.DataFrame:
    dw = dw.copy()
    dw["kind"] = np.where(dw.duration_s >= LONG_S, "long", "brief")
    # brief dwells lying between two long stalls of the same plateau, in time,
    # with the fork moving forward from the first long stall to the second
    dw["between_long"] = False
    dw["prev_long_bp"] = np.nan
    dw["next_long_bp"] = np.nan
    for _, g in dw.groupby(["data_id", "bead", "step"]):
        g = g.sort_values("k")
        longs = g[g.kind == "long"]
        for (i1, a), (i2, b) in zip(longs.iterrows(), list(longs.iterrows())[1:]):
            if b.level_bp <= a.level_bp + 30:
                continue
            mid = g[(g.k > a.k) & (g.k < b.k) & (g.kind == "brief")]
            for j in mid.index:
                dw.loc[j, ["between_long", "prev_long_bp", "next_long_bp"]] = [True, a.level_bp, b.level_bp]
    dw["from_prev_long_bp"] = dw.level_bp - dw.prev_long_bp
    return dw


def region(x):
    for name, (lo, hi) in (("N1", N1), ("linker", (N1[1], N2[0])), ("N2", N2)):
        if lo <= x < hi:
            return name
    return "before N1" if x < N1[0] else "after N2"


def shade(ax, horizontal=False):
    for lo, hi in (N1, N2):
        (ax.axhspan if horizontal else ax.axvspan)(lo, hi, color=PALETTE["light_grey"], alpha=0.35, lw=0, zorder=0)
    for x in TURN_EDGES + DYADS:
        (ax.axhline if horizontal else ax.axvline)(x, color=PALETTE["dark_grey"],
                                                   ls="-" if x in DYADS else ":", lw=0.6, zorder=0)


def figures(dw: pd.DataFrame, traces: dict):
    bins = np.arange(-100, 1101, 10)
    fig, axes = plt.subplots(3, 1, figsize=(7.0, 6.2), sharex=True)
    for ax, sel, lab, c in ((axes[0], dw.kind == "long", "long stalls (>= 1 s)", "black"),
                            (axes[1], dw.kind == "brief", "all brief dwells (0.1-1 s)", PALETTE["blue"]),
                            (axes[2], dw.between_long, "brief dwells between two long stalls", PALETTE["red"])):
        shade(ax)
        ax.hist(dw[sel].level_bp, bins=bins, color=c, alpha=0.85, lw=0)
        ax.text(0.99, 0.9, f"{lab}, n = {int(sel.sum())}", transform=ax.transAxes, ha="right", fontsize=7)
        style_axes(ax, ylabel="Count")
    axes[-1].set_xlabel("Dwell position Luz, dWLC (bp)")
    fig.savefig(OUT / "fig_dwell_positions.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_dwell_positions.png", dpi=250, bbox_inches="tight")
    plt.close(fig)

    b = dw[dw.between_long]
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.hist(b.from_prev_long_bp, bins=np.arange(0, 401, 10), color=PALETTE["red"], alpha=0.85, lw=0)
    for x in (40, 80, 120):
        ax.axvline(x, color=PALETTE["dark_grey"], ls=":", lw=0.6)
    style_axes(ax, xlabel="Brief dwell - preceding long stall (bp)", ylabel="Count")
    fig.savefig(OUT / "fig_brief_relative_to_previous_stall.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fig_brief_relative_to_previous_stall.png", dpi=250, bbox_inches="tight")
    plt.close(fig)


def example_pages(dw: pd.DataFrame, traces: dict, n_max: int = 60):
    """Plateaus with brief dwells between long stalls, largest number first."""
    key = dw[dw.between_long].groupby(["data_id", "bead", "step"]).size().sort_values(ascending=False)
    with PdfPages(OUT / "examples_brief_between_long.pdf") as pdf:
        for (d, b, s), n in key.head(n_max).items():
            tt, luz, segs = traces[(d, b, s)]
            g = dw[(dw.data_id == d) & (dw.bead == b) & (dw.step == s)]
            fig, (ax, axh) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [3, 1]}, sharey=True)
            shade(ax, horizontal=True); shade(axh, horizontal=True)
            ax.plot(tt, luz, ".", ms=1.5, color=PALETTE["dark_grey"])
            for r in g.itertuples():
                c = "black" if r.kind == "long" else (PALETTE["red"] if r.between_long else PALETTE["blue"])
                ax.plot([r.t_start_s, r.t_start_s + r.duration_s], [r.level_bp] * 2, "-", color=c, lw=2)
            axh.hist(luz, bins=np.arange(-100, 1101, 8), orientation="horizontal", color=PALETTE["dark_grey"])
            ax.set_ylim(min(-50, np.nanpercentile(luz, 1) - 30), max(700, np.nanpercentile(luz, 99) + 30))
            style_axes(ax, xlabel="Time in plateau (s)", ylabel="Luz, dWLC (bp)")
            style_axes(axh, xlabel="Frames")
            fig.suptitle(f"data_{d:03d} {b} step {s} ({g.force_pn.iloc[0]:.1f} pN): black long >= {LONG_S:.0f} s, "
                         f"red brief between long stalls, blue other brief; dotted dyad -+40, solid dyad",
                         fontsize=7.5)
            pdf.savefig(fig, dpi=200, bbox_inches="tight"); plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", type=Path, default=None)
    args = ap.parse_args()
    p.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    picks = tethers(args.selection)
    titr, _ = p.titration_files()
    frames, traces = [], {}
    for d in titr:
        dw, tr = plateau_dwells(d, picks)
        frames.append(dw); traces.update(tr)
        print(f"data_{d:03d}: {dw.groupby(['bead']).ngroups if len(dw) else 0} tethers, {len(dw)} dwells", flush=True)
    dw = classify(pd.concat(frames, ignore_index=True))
    dw["region"] = dw.level_bp.map(region)
    dw.to_csv(OUT / "dwells_all.csv", index=False)
    figures(dw, traces)
    example_pages(dw, traces)
    print(f"\n{dw.groupby(['data_id','bead']).ngroups} tethers, {dw.groupby(['data_id','bead','step']).ngroups} plateaus, "
          f"{len(dw)} dwells: {int((dw.kind=='long').sum())} long, {int((dw.kind=='brief').sum())} brief, "
          f"{int(dw.between_long.sum())} brief between two long stalls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
