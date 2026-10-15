#!/usr/bin/env python3
"""20261009 ladder: Luz model B peak map for EACH force, the user's 96 selected tethers.

User, 2026-10-11: "我希望你能分别统计每一个力的luz peakmap".

Samples: the 10 s plateau of every step the bead survived (ladder_norm1028_20261009
rules: zero = max extension of the 10 pN step, drift through the rests, steps after
a loss / fly-away / an unclosed rest dropped), from 0.3 s after the 2 mm/s move
(settling) to the end of the hold or the loss. Luz_B = (ext - zero - 1800 bp handle
stretch from 10 pN) / 2 x_ss(F), absolute bp. Each tether counts once per force
(its histogram is normalised to 1 before pooling), so long-tracked tethers do not
dominate. Histogram, 5 bp bins, square-root of the pooled density, the
basic-unzipping-analysis adaptive local-background peak calling (160 bp window,
20 % quantile, min distance 80 bp, <= 8 peaks, prominence >= 1.5 robust sigma).

Outputs (WHZ_20261009/luz_modelB/):
  fig_luzB_peakmap_by_force.(pdf|png|svg)   one row per force + a force x Luz map
  luzB_peaks_by_force.csv                    every peak: force, bp, height, nearest 601 landmark
  luzB_density_by_force.csv                  the pooled density per force (for re-plotting)

    .venv/bin/python luzB_peakmap_by_force_20261009.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MODEL", "B")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

import ladder_norm1028_20261009 as L  # noqa: E402
import process_20261009 as p  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402

# SCALE=norm1028 (user, 2026-10-11: the 0-1028 trace "更正确"): the user's per-tether
# normalisation 1028 * (ext - zero) / span instead of the absolute Luz model B
SCALE = os.environ.get("SCALE", "luzB")
OUT = p.DEST / ("norm1028_by_force" if SCALE == "norm1028" else "luz_modelB")
TAG = "norm1028" if SCALE == "norm1028" else "luzB"
XLAB = ("Normalised position (0-1028, per tether)" if SCALE == "norm1028"
        else "Luz, model B (bp from the 10 pN zero)")
SETTLE_S = 0.3
EDGES = np.arange(-100.0, 1150.0 + 5.0, 5.0)
LANDMARKS = {"N1 entry": 314, "N1 dyad-40": 347, "N1 dyad": 387, "N1 dyad+40": 427, "N1 exit": 460,
             "N2 entry": 511, "N2 dyad-40": 544, "N2 dyad": 584, "N2 dyad+40": 624, "N2 exit": 657,
             "trunk end": 971}


def collect() -> dict[float, list[np.ndarray]]:
    sel = pd.read_csv(L.SEL_CSV)
    per_force: dict[float, list[np.ndarray]] = {}
    for d, g in sel.groupby("data_id"):
        run = p.load(int(d))
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{int(d):03d}.h5")
        t = run["time_s"]
        for bead in g.bead:
            roi = bead[1:]
            ext, rec = L.bead_scale(run, steps, roi, trk.get(roi, {}))
            if not np.isfinite(rec["zero_nm"]):
                continue
            for st in steps:
                if st["step"] == 1 or st["step"] not in rec["usable_steps"]:
                    continue
                pl = st["plateau"]
                keep = (t[pl] - t[pl[0]] >= SETTLE_S) & (t[pl] < rec["t_loss"]) & np.isfinite(ext[pl])
                if keep.sum() < 20:
                    continue
                f = np.full(keep.sum(), st["force_pn"])
                if SCALE == "norm1028":
                    if not (np.isfinite(rec["span_nm"]) and rec["span_nm"] > 0):
                        continue
                    luz = 1028.0 * (ext[pl][keep] - rec["zero_nm"]) / rec["span_nm"]
                else:
                    luz = L.model_luz(ext[pl][keep] - rec["zero_nm"], f)
                per_force.setdefault(round(st["force_pn"], 1), []).append(luz[np.isfinite(luz)])
        print(f"data_{int(d):03d}", flush=True)
    return per_force


def peak_map(samples: list[np.ndarray]):
    a = L.peak_args()
    centers = 0.5 * (EDGES[:-1] + EDGES[1:])
    h = np.zeros(centers.size)
    for x in samples:
        c, _ = np.histogram(x, bins=EDGES)
        if c.sum():
            h += c / c.sum()
    sq = np.sqrt(h)
    smax = sq.max() if sq.max() > 0 else 1.0
    sm = L.BUA.centered_mean(sq, a.hist_smooth_bins)
    base, sig = L.BUA.local_peak_signal(centers, sm, a.peak_local_window_bp, a.peak_baseline_quantile)
    peaks = L.BUA.find_histogram_peaks(centers, sm, sig, sm / smax, "adaptive", a.peak_min_height,
                                       a.peak_min_prominence * smax, a.peak_noise_sigma,
                                       a.peak_min_distance_bp, -100.0, 8)
    return centers, sq / smax, base / smax, sig / smax, peaks, h


def nearest(x: float) -> str:
    if SCALE == "norm1028":
        LANDMARKS["trunk end"] = 1028
    k = min(LANDMARKS, key=lambda n: abs(LANDMARKS[n] - x))
    return f"{k} ({x - LANDMARKS[k]:+.0f})"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p.apply_style()
    pf = collect()
    forces = [f for f in sorted(pf) if f >= 18.0 and len(pf[f]) >= 5]
    rows, dens = [], {}
    fig = plt.figure(figsize=(7.09, 1.05 * len(forces) + 2.6))
    gs = fig.add_gridspec(len(forces) + 1, 1, height_ratios=[1.0] * len(forces) + [2.2], hspace=0.08, top=0.975, bottom=0.05)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(len(forces))]
    mapax = fig.add_subplot(gs[-1, 0], sharex=axes[0])
    img = []
    for ax, f in zip(axes, forces):
        c, sq, base, sig, peaks, h = peak_map(pf[f])
        dens[f] = h
        img.append(sq)
        for lo, hi in p.NUCLEOSOME_SITES_BP:
            ax.axvspan(lo, hi, color=PALETTE["light_grey"], alpha=0.5, lw=0, zorder=0)
        for x0 in (387, 584):
            ax.axvline(x0, color=PALETTE["dark_grey"], lw=0.5, ls=":", zorder=1)
        ax.axvline(1028 if SCALE == "norm1028" else 971, color="black", lw=0.6, ls="--")
        ax.bar(c, sq, width=5.0, color=PALETTE["grey"], lw=0, alpha=0.55)
        ax.plot(c, base, color=PALETTE["dark_grey"], lw=0.6, ls="--")
        ax.plot(c, sig, color=PALETTE["blue"], lw=0.9)
        for x in peaks:
            ax.plot(x, 1.0, "v", ms=3.5, color=PALETTE["red"], clip_on=False, zorder=5)
            ax.text(x, 1.04, f"{x:.0f}", color=PALETTE["red"], fontsize=5.5, ha="center", va="bottom")
            rows.append({"force_pn": f, "peak_bp_B": round(x, 0), "scale": TAG, "n_tethers": len(pf[f]),
                         "signal": round(float(np.interp(x, c, sig)), 3), "nearest_landmark": nearest(x)})
        ax.text(1.005, 0.5, f"{f:.1f} pN\nn = {len(pf[f])}", transform=ax.transAxes, fontsize=6.5,
                va="center", ha="left")
        ax.set_ylim(0, 1.18)
        ax.set_yticks([0, 1])
        style_axes(ax)
        ax.tick_params(labelbottom=False, labelsize=6)
    im = np.array(img)
    c = 0.5 * (EDGES[:-1] + EDGES[1:])
    mapax.imshow(im, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest",
                 extent=(EDGES[0], EDGES[-1], -0.5, len(forces) - 0.5))
    for r in rows:
        mapax.plot(r["peak_bp_B"], forces.index(r["force_pn"]), "o", ms=2.5, mfc=PALETTE["red"], mec="none")
    for x0 in (314, 460, 511, 657):
        mapax.axvline(x0, color=PALETTE["blue"], lw=0.6, ls="--")
    mapax.set_yticks(range(len(forces)), [f"{f:.1f}" for f in forces], fontsize=6)
    style_axes(mapax, xlabel=XLAB, ylabel="Force (pN)")
    mapax.set_xlim(EDGES[0], EDGES[-1])
    fig.text(0.04, 0.62, "sqrt pooled density (norm.)", rotation=90, va="center", fontsize=8)
    fig.suptitle(f"20261009, your 96 selected tethers: {TAG} peak map per force (plateaus, from "
                 f"{SETTLE_S:.1f} s; each tether weight 1; 601 sites shaded, dyads dotted, " + ("1028" if SCALE == "norm1028" else "trunk end 971") + " dashed; "
                 "red = peaks)", fontsize=7.5, y=0.995)
    for ext_ in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"fig_{TAG}_peakmap_by_force.{ext_}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pk = pd.DataFrame(rows)
    pk.to_csv(OUT / f"{TAG}_peaks_by_force.csv", index=False)
    pd.DataFrame({"bp_center": c, **{f"{f:.1f}pN": dens[f] for f in forces}}).to_csv(
        OUT / f"{TAG}_density_by_force.csv", index=False)
    (OUT / f"fig_{TAG}_peakmap_by_force.txt").write_text(
        "Luz model B peak map per force, your 96 selected tethers (2 x_ss per opened bp, zero = max "
        "extension of the 10 pN step). Each row: the plateau samples of every step at that force that the "
        "bead survived, from 0.3 s after the move; each tether normalised to weight 1; 5 bp bins; grey bars "
        "= square root of the pooled density (normalised to its maximum); dashed = local 20 % quantile "
        "background (160 bp window); blue = background-corrected signal; red triangles and numbers = "
        "adaptive peaks (min distance 80 bp, prominence >= 1.5 robust sigma). Grey bands: 601 sites "
        "314-460 / 511-657 bp; dotted: dyads 387 / 584; dashed: trunk end 971. Bottom: the same densities "
        "as a force x Luz map, red dots = peaks, blue dashed = 601 site edges.\n")
    print(pk.to_string(index=False))
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
