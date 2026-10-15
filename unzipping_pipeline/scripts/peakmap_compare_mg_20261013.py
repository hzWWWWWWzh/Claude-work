#!/usr/bin/env python3
"""Unzipping stalling peak maps per force: -Mg2+ (20261009, 96 tethers) vs +1 mM Mg2+ (20261013, 101 tethers).

User (2026-10-14): "给出-Mg/+Mg的Unzipping stalling peak 对比图" (with the 20261009 Luz model B per-force
peak map as the template).

Same method as luzB_peakmap_by_force_20261009.py, reused unchanged (collect / peak_map): samples = the
plateau of every step the bead survived, from 0.3 s after the move; Luz_B = (ext - zero - 1800 bp handle
stretch from 10 pN) / 2 x_ss(F) (absolute bp from the 10 pN zero; independent of the 1028 point, so the
46 pN reference step of 20261013 does not enter); each tether weight 1; 5 bp bins; sqrt of the pooled
density; adaptive local-background peaks (basic-unzipping-analysis settings). The two batches are collected
by re-pointing process_20261009.RAW / DEST and the selection list between the two runs.

Outputs (Proecessed Traces/WHZ_20261013/compare_vs_1009/peakmap_mg/):
  fig_peakmap_mg_by_force.(pdf|png|svg)   one row per force, -Mg blue / +Mg red outlines + peaks
  fig_peakmap_mg_heat.(pdf|png|svg)       force x Luz_B heat maps side by side
  peaks_mg_by_force.csv, density_mg_by_force.csv

    .venv/bin/python peakmap_compare_mg_20261013.py
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

os.environ["MODEL"] = "B"
os.environ["SCALE"] = "luzB"
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import process_20261009 as p  # noqa: E402
import ladder_norm1028_20261009 as L  # noqa: E402
import luzB_peakmap_by_force_20261009 as M  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402

ROOT = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces")
BATCHES = {
    "20261009": {"label": "-Mg2+ (20261009)", "sel": ROOT / "WHZ_20261009/selection_96_tethers.csv",
                 "colour": PALETTE["blue"]},
    "20261013": {"label": "+1 mM Mg2+ (20261013)", "sel": ROOT / "WHZ_20261013/selection_101_tethers.csv",
                 "colour": PALETTE["red"]},
}
OUT = ROOT / "WHZ_20261013/compare_vs_1009/peakmap_mg"
CACHE = OUT / "_samples.pkl"
FORCES = (18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9)


def collect_batch(batch: str) -> dict:
    p.RAW = Path(f"/Users/whz/Desktop/DATA MT/unzipping/{batch}")
    p.DEST = ROOT / f"WHZ_{batch}"
    L.SEL_CSV = BATCHES[batch]["sel"]
    return M.collect()


def load(recollect=False) -> dict:
    if CACHE.exists() and not recollect:
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)
    OUT.mkdir(parents=True, exist_ok=True)
    res = {b: collect_batch(b) for b in BATCHES}
    with open(CACHE, "wb") as fh:
        pickle.dump(res, fh)
    return res


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p.apply_style()
    data = load("--recollect" in sys.argv)
    centers = 0.5 * (M.EDGES[:-1] + M.EDGES[1:])
    rows, dens = [], {}
    fig = plt.figure(figsize=(7.09, 1.0 * len(FORCES) + 0.9))
    gs = fig.add_gridspec(len(FORCES), 1, hspace=0.1, top=0.95, bottom=0.07)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(len(FORCES))]
    for ax, f in zip(axes, FORCES):
        for lo, hi in p.NUCLEOSOME_SITES_BP:
            ax.axvspan(lo, hi, color=PALETTE["light_grey"], alpha=0.5, lw=0, zorder=0)
        for x0 in (387, 584):
            ax.axvline(x0, color=PALETTE["dark_grey"], lw=0.5, ls=":", zorder=1)
        ax.axvline(971, color="black", lw=0.6, ls="--")
        txt = [f"{f:g} pN"]
        for k, (b, info) in enumerate(BATCHES.items()):
            s = data[b].get(f, [])
            if len(s) < 5:
                txt.append(f"n{'-' if k == 0 else '+'} = {len(s)}")
                continue
            c, sq, base, sig, peaks, h = M.peak_map(s)
            dens[(b, f)] = h
            ax.step(c, sq, where="mid", color=info["colour"], lw=0.8, alpha=0.9)
            for x in peaks:
                yv = 1.08 if k == 0 else 1.0
                ax.plot(x, yv, "v", ms=3.2, color=info["colour"], clip_on=False, zorder=5)
                rows.append({"batch": b, "force_pn": f, "peak_luzB_bp": round(x, 0), "n_tethers": len(s),
                             "signal": round(float(np.interp(x, c, sig)), 3)})
            txt.append(f"n{'-' if k == 0 else '+'} = {len(s)}")
        ax.text(1.005, 0.5, "\n".join(txt), transform=ax.transAxes, fontsize=6.5, va="center", ha="left")
        ax.set_ylim(0, 1.15)
        ax.set_yticks([0, 1])
        ax.set_xlim(M.EDGES[0], M.EDGES[-1])
        style_axes(ax)
        ax.tick_params(labelbottom=False, labelsize=6)
    axes[-1].tick_params(labelbottom=True)
    axes[-1].set_xlabel("Luz, model B (bp from the 10 pN zero)", fontsize=8)
    fig.text(0.035, 0.5, "sqrt pooled density (norm.)", rotation=90, va="center", fontsize=8)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color=i["colour"], lw=1.0, label=i["label"]) for i in BATCHES.values()],
               loc="upper center", ncol=2, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.995))
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"fig_peakmap_mg_by_force.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # heat maps, same colour scale
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.6), sharey=True)
    for ax, (b, info) in zip(axes, BATCHES.items()):
        img = np.array([np.sqrt(dens[(b, f)]) / np.sqrt(dens[(b, f)]).max() if (b, f) in dens else np.full(centers.size, np.nan)
                        for f in FORCES])
        ax.imshow(img, aspect="auto", origin="lower", cmap="Greys", interpolation="nearest",
                  extent=(M.EDGES[0], M.EDGES[-1], -0.5, len(FORCES) - 0.5))
        for r in [r for r in rows if r["batch"] == b]:
            ax.plot(r["peak_luzB_bp"], FORCES.index(r["force_pn"]), "o", ms=2.2, mfc=info["colour"], mec="none")
        for x0 in (314, 460, 511, 657):
            ax.axvline(x0, color=PALETTE["blue"], lw=0.5, ls="--")
        ax.set_yticks(range(len(FORCES)), [f"{f:g}" for f in FORCES], fontsize=6)
        style_axes(ax, xlabel="Luz, model B (bp)", ylabel="Force (pN)" if b == "20261009" else None)
        ax.set_title(info["label"], fontsize=7, loc="left")
        ax.tick_params(labelsize=6)
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"fig_peakmap_mg_heat.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    pk = pd.DataFrame(rows)
    pk.to_csv(OUT / "peaks_mg_by_force.csv", index=False)
    pd.DataFrame({"luzB_bp_center": centers, **{f"{b}_{f:g}pN": dens[(b, f)] for (b, f) in dens}}).to_csv(
        OUT / "density_mg_by_force.csv", index=False)
    print(pk.pivot_table(index="force_pn", columns="batch", values="peak_luzB_bp",
                         aggfunc=lambda v: "; ".join(f"{x:.0f}" for x in sorted(v))).to_string())
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
