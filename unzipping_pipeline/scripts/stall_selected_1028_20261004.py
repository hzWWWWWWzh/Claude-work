#!/usr/bin/env python3
"""20261004: stall statistics on the tethers the user selected.

Input: the user's current selection (66 tethers left in the trimmed merged PDF,
column "selected" of norm1028_v6/norm1028_beads.csv) and the v6 normalisation (0 = each file's C1 extension at 10 pN; v5.3 rule below)
(normalise_1028_20261004.py): zero = the closed molecule fully stretched, one
closed stretch per tether on top of each file's own 5 pN level (lower of
C1/C2, stuck cycles left out); a file whose ramp re-zips started open and takes
the 5 pN level of the nearest earlier re-closed file (ZERO_MODE=own: keeps its
own, the sensitivity check); 1028 = the furthest the fork got on the complete Fmax
plateaus of the highest-force file with a complete pull.

Stall at a given Fmax = median over the last 1 s of that cycle's ~10 s plateau.
A cycle counts ("used") if the hold lasted >= 9 s (FOV02's 1-3 s holds are not
comparable), the bead is tracked to the end of the hold, the file is not a LUT
level jump, and the fork is stationary in that last second: IQR of the
extension <= 75 nm (about 100 bp at 30 pN; judged in nm because a small span
multiplies noise on the 1028 scale). Per tether and force: median of the used
cycles.

Two scales:
  * 1028 scale, the user's rule (v5 zero and end, above). Each tether's own
    1028 file is left out of the statistics (its maximum, not its plateau
    end, is fixed); tethers with a span < 150 nm are left out of this scale
    (x > 7 noise; verify6: FOV07 R040 alone moved the 39.9 pN LS shift from
    481 to 119). Tethers whose closed state is seen only below 10 pN (no
    usable 0 point, v5.1) are left out of both scales. Rule B (largest maximum over complete pulls) is kept as
    columns in stall_by_cycle.csv only.
  * dWLC bp, span-free: bp = (extension above the zero - handle stretch from
    the zero force F* to F) / (nm per opened bp at F), with the same zeros. The
    physical estimate of how far the fork got; it needs no 1028 point, so
    every file is used. On the 1028 scale, "1028" is each tether's own
    maximum, a full 1028 bp unzip only for full-length tethers (span within
    0.80-1.25 of the dWLC full-unzip span).

Statistics (ladders are incomplete, so all available pairs are used): a
median-type (least absolute deviations) additive model stall = tether + force
level (levels with >= 5 tethers; anchor = median tether intercept; the LAD
optimum is not unique, shifts are indeterminate by ~25 on the 1028 scale),
95% CI by resampling tethers (one bootstrap, reused everywhere); neighbouring
forces of each field's ladder on tethers that have both (Wilcoxon,
Bonferroni); per tether highest vs lowest force and Spearman rho; C2 vs C1 in
one file; 601 sites on full-length tethers with the tether as the unit
(Wilcoxon and sign test, windows shifted by -100..+200 bp, dWLC coordinates).

Caveats built in:
  * the user selected these tethers by eye from pages that showed the whole
    titration and the 601 bands, so the statistics describe the selected
    set; the unselected normalised tethers of the same fields are analysed
    the same way for comparison (selection_comparison.csv);
  * Fmax rises with acquisition time in every field, so force and pulling
    history cannot be separated.

Outputs: WHZ_20261004/norm1028_v5_3[_ownbase|_beadbase]/stall_selected/ (each figure with a caption .txt)
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
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, NullLocator
from scipy.optimize import linprog
from scipy.stats import binomtest, spearmanr, wilcoxon
from whz_plot import PALETTE, apply_style, export_panel, raw_trace, style_axes

from normalise_1028_20261004 import (
    MIN_SPAN_STATS_NM, NORM_BP, OUT as V2, bead_scale, complete_plateaus, extension_nm,
    file_status, last_plateau_end, plateau_max,
)
from process_20261004 import CALIBRATION, DEST, FOVS, NUCLEOSOME_SITES_BP, RAW
from reprocess_26_core_sop import HANDLE_BP, delta_nm_per_bp, x_ds
from tether_stall_titration_20261004 import plateaus
from unzip_raw_points import load_raw_run

OUT = V2 / "stall_selected"
REQUIRED_HOLD_S = 9.0
END_WINDOW_S = 1.0
END_REACH_S = 0.1
MAX_END_IQR_NM = 75.0
ZERO_PN = 15.0
FULL_RATIO = (0.80, 1.25)
MIN_TETHERS_PER_LEVEL = 5
N_BOOT = 1000
YLIM = (-300.0, 1260.0)          # one limit for every figure (30-bin edges end here)
LEVELS_ALL = [18.0, 25.9, 27.9, 30.0, 34.6, 39.9, 42.8, 46.0]
SHIFTS_BP = (-100, -50, -25, 0, 25, 50, 100, 200)
# The user's final rule (2026-10-08): rule A with a robust zero - the 1028 file
# uses the bead's median zero when its own is > 150 nm off (bead_scale in
# normalise_1028_20261004.py). Rule B columns stay in stall_by_cycle.csv for the
# record but are no longer analysed.
SCALES = {
    "A": ("norm1028_scale", "stall_A", "is_1028_file_A", "span_ratio_A",
          "Stall (1028 scale)"),
    "bp": ("dwlc_bp", "stall_bp", None, "span_ratio_A", "Stall (bp, dWLC)"),
}


def full_unzip_span_nm(f_pn: float, zero_pn: float = ZERO_PN) -> float:
    """dWLC extension gain of a fully unzipped 1028 bp construct, zero -> f."""
    f = np.array([f_pn])
    return float(HANDLE_BP * (x_ds(f) - x_ds(np.array([zero_pn])))[0]
                 + NORM_BP * delta_nm_per_bp(f)[0])


def dwlc_bp(nm: float, f_pn: float, f_zero_pn: float = ZERO_PN) -> float:
    """Opened bp from the extension above the zero, at force f. The zero is the
    closed molecule at f_zero (v5: F*, the highest force still closed), so the
    handle's WLC stretch from f_zero to f is taken off first."""
    f = np.array([f_pn])
    handle = HANDLE_BP * (x_ds(f) - x_ds(np.array([f_zero_pn])))[0]
    return float((nm - handle) / delta_nm_per_bp(f)[0])


def complete_files(files, runs, stats, rec):
    """Files with a complete pull (same rule as bead_scale), up to the break."""
    return list(rec.get("complete_files", []))


def collect(sel: pd.DataFrame, first: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fov, beads in sel.groupby("fov").bead:
        files = FOVS[fov]
        runs = {d: load_raw_run(RAW, d, CALIBRATION) for d in files}
        fa = first[first.fov == fov]
        for bead in beads:
            roi = bead[1:]
            stats, exts = {}, {}
            for d in files:
                r = fa[(fa.data_id == d) & (fa.bead == bead)]
                if r.empty:
                    continue
                stats[d] = file_status(r.iloc[0])
                if stats[d]["present"] and roi in runs[d].z_um:
                    tl = stats[d]["t_loss"]
                    exts[d], _ = extension_nm(runs[d].z_um[roi], runs[d].time_s,
                                              runs[d].force_pn, tl if np.isfinite(tl) else None)
            rec = bead_scale(files, runs, exts, stats)
            if not rec["zeros"]:
                continue
            if rec.get("zero_unreliable"):
                # closed only below 10 pN: no usable 0 point (v5.1)
                print(f"  {fov} {bead}: left out, 0 point unusable (F* {rec['f_star']:.0f} pN)")
                continue
            # a bead without a usable span (e.g. span <= 0 after the zero rule)
            # still has dWLC bp values; only its 1028-scale values are NaN
            cand = {}
            for d in complete_files(files, runs, stats, rec):
                zm, _, _ = plateau_max(runs[d], exts[d], stats[d],
                                       complete_plateaus(runs[d], stats[d]))
                if np.isfinite(zm):
                    cand[d] = zm - rec["zeros"][d][0]
            ref_a = int(rec["max_file"]) if np.isfinite(rec["max_file"]) else -1
            ref_b = max(cand, key=cand.get) if cand else ref_a
            span_b = cand.get(ref_b, rec["span_nm"])
            # the 1028 scale needs a span of >= 150 nm (else noise x > 7)
            spans = {"A": (rec["span_nm"] if rec["normalised"]
                           and rec["span_nm"] >= MIN_SPAN_STATS_NM else np.nan, ref_a),
                     "B": (span_b, int(ref_b))}
            f_ref = {k: round(float(runs[v[1]].force_pn.max()), 1) if v[1] in runs else np.nan
                     for k, v in spans.items()}
            for d, z in rec["zeros"].items():
                st, run = stats[d], runs[d]
                e_nm = exts[d] - z[0]                     # extension above this file's zero
                fmax = round(float(run.force_pn.max()), 1)
                for c, seg in plateaus(run).items():
                    t = run.time_s[seg]
                    hold = float(t[-1] - t[0])
                    keep = np.isfinite(e_nm[seg])
                    if np.isfinite(st["t_loss"]):
                        keep &= t < st["t_loss"]
                    if keep.sum() < 10:
                        continue
                    tr, v = t[keep] - t[0], e_nm[seg][keep]
                    reached = tr[-1] >= hold - END_REACH_S
                    end = v[tr >= hold - END_WINDOW_S] if reached else v[:0]
                    if end.size >= 3:
                        q1, q3 = np.percentile(end, [25, 75])
                        val = float(np.median(end))
                    else:
                        q1 = q3 = val = np.nan
                    last3 = tr >= hold - 3.0
                    slope = (float(np.polyfit(tr[last3], v[last3], 1)[0])
                             if reached and last3.sum() >= 10 else np.nan)
                    why = ("hold < 9 s" if hold < REQUIRED_HOLD_S else
                           "level-jump file" if st["jump"] else
                           "lost before the end of the hold" if not reached else
                           "not stationary (IQR > 75 nm)" if not (q3 - q1) <= MAX_END_IQR_NM else "")
                    row = {"fov": fov, "bead": bead, "tether": f"{fov}_{bead}", "data_id": d,
                           "fmax_pn": fmax, "cycle": c, "hold_s": round(hold, 2),
                           "stall_nm": round(val, 1),
                           "stall_bp": round(dwlc_bp(val, fmax, z[5]), 1) if np.isfinite(val) else np.nan,
                           "zero_force_pn": round(z[5], 1),
                           "end_iqr_nm": round(float(q3 - q1), 1),
                           "end_slope_nm_per_s": round(slope, 1),
                           "used": why == "", "why_not_used": why}
                    for k, (span, ref) in spans.items():
                        ok_span = np.isfinite(span) and span > 0
                        row[f"stall_{k}"] = round(NORM_BP * val / span, 1) if ok_span else np.nan
                        row[f"span_{k}_nm"] = round(span, 1) if np.isfinite(span) else np.nan
                        row[f"ref_file_{k}"] = ref
                        row[f"is_1028_file_{k}"] = d == ref
                        row[f"span_ratio_{k}"] = (round(span / full_unzip_span_nm(f_ref[k]), 3)
                                                  if ok_span and np.isfinite(f_ref[k]) else np.nan)
                    rows.append(row)
        print(f"{fov}: {beads.size} tethers", flush=True)
    return pd.DataFrame(rows)


def scale_view(st: pd.DataFrame, key: str) -> pd.DataFrame:
    """stall / span_ratio / is_1028_file columns for one scale; on the dWLC bp
    scale no file is pinned, so none is left out."""
    _, val, ref, ratio, _ = SCALES[key]
    v = st.assign(stall=st[val], span_ratio=st[ratio])
    v["is_1028_file"] = st[ref] if ref else False
    return v[v.stall.notna() | ~v.used]


def levels_table(st: pd.DataFrame) -> pd.DataFrame:
    """One value per tether and force: median of the used cycles that are not
    in the tether's 1028 file (files sharing an Fmax are never merged with it)."""
    u = st[st.used & ~st.is_1028_file]
    return (u.groupby(["tether", "fov", "bead", "fmax_pn"], as_index=False)
             .agg(stall=("stall", "median"), n_cycles=("cycle", "nunique"),
                  span_ratio=("span_ratio", "first")))


def ref_levels(st: pd.DataFrame) -> pd.DataFrame:
    """The 1028 file's own plateau-end values, for display only."""
    u = st[st.used & st.is_1028_file]
    return u.groupby(["tether", "fmax_pn"], as_index=False).agg(stall=("stall", "median"))


def design(t: pd.DataFrame, levels: list[float]):
    tethers = sorted(t.tether.unique())
    ti = {k: i for i, k in enumerate(tethers)}
    li = {f: j for j, f in enumerate(levels)}
    x = np.zeros((len(t), len(tethers) + len(levels) - 1))
    for r, (k, f) in enumerate(zip(t.tether, t.fmax_pn)):
        x[r, ti[k]] = 1.0
        if li[f]:
            x[r, len(tethers) + li[f] - 1] = 1.0
    return x, len(tethers)


def fit_shift(t: pd.DataFrame, levels: list[float], lad: bool = True):
    """Level effects relative to the lowest level within tethers, and the
    MEDIAN tether intercept as the anchor (robust to one runaway tether)."""
    present = [f for f in levels if f in set(t.fmax_pn)]
    x, nt = design(t, present)
    y = t.stall.to_numpy(float)
    if not lad:
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
    else:
        n, p = x.shape
        c = np.r_[np.zeros(p), np.ones(2 * n)]
        a_eq = np.hstack([x, np.eye(n), -np.eye(n)])
        bounds = [(None, None)] * p + [(0, None)] * (2 * n)
        beta = linprog(c, A_eq=a_eq, b_eq=y, bounds=bounds, method="highs").x[:p]
    eff = np.r_[0.0, beta[nt:]]
    out = {f: float(eff[j]) for j, f in enumerate(present)}
    return {f: out.get(f, np.nan) for f in levels}, float(np.median(beta[:nt]))


def boot_shift(t: pd.DataFrame, levels: list[float], lad: bool = True, n_boot: int = N_BOOT):
    rng = np.random.default_rng(0)
    groups = {k: g for k, g in t.groupby("tether")}
    keys = list(groups)
    res = []
    for _ in range(n_boot):
        pick = rng.choice(len(keys), len(keys), replace=True)
        bt = pd.concat([groups[keys[i]].assign(tether=f"{keys[i]}#{j}")
                        for j, i in enumerate(pick)], ignore_index=True)
        if bt.fmax_pn.min() != levels[0]:
            continue
        res.append(fit_shift(bt, levels, lad)[0])
    return pd.DataFrame(res)


def model_set(t: pd.DataFrame):
    n_lv = t.groupby("fmax_pn").tether.nunique()
    levels = sorted(n_lv[n_lv >= MIN_TETHERS_PER_LEVEL].index)
    tm = t[t.fmax_pn.isin(levels)]
    tm = tm[tm.groupby("tether").fmax_pn.transform("nunique") >= 2]
    return tm, levels


def step_tests(t: pd.DataFrame) -> pd.DataFrame:
    rows = {}
    for (fov, tether), g in t.groupby(["fov", "tether"]):
        lad = sorted(set(t[t.fov == fov].fmax_pn))
        s = g.set_index("fmax_pn").stall
        for a, b in zip(lad[:-1], lad[1:]):
            if a in s.index and b in s.index:
                rows.setdefault((a, b), []).append(s[b] - s[a])
    out = []
    for (a, b), d in sorted(rows.items()):
        d = np.array(d)
        p = wilcoxon(d).pvalue if len(d) >= 5 and (d != 0).any() else np.nan
        out.append({"comparison": f"{a}->{b} pN", "n": len(d), "median": round(np.median(d), 0),
                    "q25": round(np.percentile(d, 25), 0), "q75": round(np.percentile(d, 75), 0),
                    "frac_up": round((d > 0).mean(), 2), "wilcoxon_p": p})
    o = pd.DataFrame(out)
    o["p_bonferroni"] = (o.wilcoxon_p * o.wilcoxon_p.notna().sum()).clip(upper=1)
    return o


def trend_tests(t: pd.DataFrame) -> pd.DataFrame:
    hl, rho = [], []
    for _, g in t.groupby("tether"):
        g = g.sort_values("fmax_pn")
        if len(g) < 2:
            continue
        hl.append(g.stall.iloc[-1] - g.stall.iloc[0])
        if len(g) >= 3:
            r = spearmanr(g.fmax_pn, g.stall)[0]
            if np.isfinite(r):
                rho.append(r)
    hl, rho = np.array(hl), np.array(rho)
    return pd.DataFrame([
        {"comparison": "highest - lowest force, per tether", "n": len(hl),
         "median": round(np.median(hl), 0), "q25": round(np.percentile(hl, 25), 0),
         "q75": round(np.percentile(hl, 75), 0), "frac_up": round((hl > 0).mean(), 2),
         "wilcoxon_p": wilcoxon(hl).pvalue},
        {"comparison": "Spearman rho stall vs force (>= 3 forces; unitless)",
         "n": len(rho), "median": round(np.median(rho), 2), "q25": round(np.percentile(rho, 25), 2),
         "q75": round(np.percentile(rho, 75), 2), "frac_up": round((rho > 0).mean(), 2),
         "wilcoxon_p": wilcoxon(rho).pvalue}])


def inside(v: float, shift: float = 0.0) -> bool:
    return any(lo + shift <= v <= hi + shift for lo, hi in NUCLEOSOME_SITES_BP)


def site_tests(t: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Full-length tethers, stalls in 0..0.9 x 1028. Unit = tether: fraction of
    its stalls inside the sites vs the fraction of the range they cover;
    Wilcoxon and sign test; the same with both sites shifted together."""
    reach = 0.9 * NORM_BP
    v = t[(t.stall >= 0) & (t.stall < reach)]
    rows = []
    for sh in SHIFTS_BP:
        cover = sum(max(0.0, min(hi + sh, reach) - max(lo + sh, 0.0))
                    for lo, hi in NUCLEOSOME_SITES_BP) / reach
        frac = v.groupby("tether").stall.apply(lambda x: x.apply(inside, shift=sh).mean())
        if len(frac) < 3:
            continue
        d = frac.values - cover
        rng = np.random.default_rng(0)
        boot = [rng.choice(frac.values, len(frac), replace=True).mean() for _ in range(5000)]
        rows.append({"coordinates": tag, "site_shift_bp": sh, "n_tethers": len(frac),
                     "n_stalls": len(v), "fraction_inside": round(frac.mean(), 2),
                     "ci_lo": round(np.percentile(boot, 2.5), 2),
                     "ci_hi": round(np.percentile(boot, 97.5), 2), "cover": round(cover, 2),
                     "excess": round(frac.mean() - cover, 2),
                     "wilcoxon_p_greater": wilcoxon(d, alternative="greater").pvalue
                     if (d != 0).any() else np.nan,
                     "sign_p_greater": binomtest(int((d > 0).sum()), int((d != 0).sum()), 0.5,
                                                 alternative="greater").pvalue
                     if (d != 0).any() else np.nan})
    return pd.DataFrame(rows)


def shade(ax, horizontal=True):
    for lo, hi in NUCLEOSOME_SITES_BP:
        (ax.axhspan if horizontal else ax.axvspan)(lo, hi, color="#EBEBEB", lw=0, zorder=0)


def fig_trend(t, refs, levels, shift, base, ylabel, out, show_sites, unit):
    """Top: every tether's stall per force (lines join one tether; open = its
    1028 file, excluded; triangles = beyond the axis). Bottom: the
    tether-adjusted shift from the lowest force with its 95% CI."""
    lv = [f for f in LEVELS_ALL if f in set(t.fmax_pn) | set(refs.fmax_pn)]
    pos = {f: k for k, f in enumerate(lv)}
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(3.39, 4.6), sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08})
    if show_sites:
        shade(ax)
    over_y = YLIM[1] + 0.03 * (YLIM[1] - YLIM[0])
    allp = pd.concat([t.assign(ref=False), refs.assign(ref=True)], ignore_index=True)
    tethers = sorted(allp.tether.unique())
    jit = dict(zip(tethers, np.random.default_rng(0).uniform(-0.15, 0.15, len(tethers))))
    ycl = lambda y: np.where(y > YLIM[1], over_y, np.maximum(y, YLIM[0]))
    for k, g in allp.groupby("tether"):
        g = g.sort_values("fmax_pn")
        ax.plot(g.fmax_pn.map(pos) + jit[k], ycl(g.stall), color="#D0D0D0", lw=0.35, alpha=0.7,
                zorder=1, rasterized=True)
    xj = allp.fmax_pn.map(pos) + allp.tether.map(jit)
    inr = allp.stall <= YLIM[1]
    for sel, filled in ((~allp.ref & inr, True), (allp.ref & inr, False)):
        raw_trace(ax, xj[sel], allp.stall[sel], color=PALETTE["dark_grey"], filled=filled,
                  marker_area=6, alpha=0.8)
    over = ~inr
    if over.any():
        raw_trace(ax, xj[over], np.full(over.sum(), over_y), marker="^", filled=True,
                  marker_area=8, color=PALETTE["dark_grey"], clip_on=False)
    s = shift.set_index("fmax_pn")
    s = s[s.index.isin(levels)]
    ax.plot([pos[f] for f in s.index], base + s.shift_lad, ls="none", marker="D", ms=4.5,
            mfc="black", mec="white", mew=0.6, zorder=6)
    ax.set_ylim(*YLIM)
    dg = PALETTE["dark_grey"]
    entries = [(Line2D([], [], ls="none", marker="o", ms=3, mfc=dg, mec=dg), "analysed")]
    if len(refs):
        entries.append((Line2D([], [], ls="none", marker="o", ms=3, mfc="none", mec=dg),
                        "1028 file (excluded)"))
    if over.any():
        entries.append((Line2D([], [], ls="none", marker="^", ms=3.5, mfc=dg, mec=dg),
                        f"> {YLIM[1]:.0f} (n = {int(over.sum())})"))
    entries.append((Line2D([], [], ls="none", marker="D", ms=4, mfc="black", mec="white"),
                    "LAD level estimate"))
    ax.legend(handles=[h for h, _ in entries], labels=[l for _, l in entries],
              loc="lower center", bbox_to_anchor=(0.5, 1.07), ncol=2, fontsize=8,
              borderaxespad=0)
    bx.axhline(0, color=PALETTE["light_grey"], lw=0.7, zorder=0)
    s1 = s[s.index != levels[0]]
    bx.errorbar([pos[f] for f in s1.index], s1.shift_lad,
                yerr=[s1.shift_lad - s1.lad_lo, s1.lad_hi - s1.shift_lad],
                fmt="D", ms=4, color="black", capsize=2, label="LAD shift, 95% CI (tethers)")
    bx.plot([pos[levels[0]]], [0], ls="none", marker="D", ms=4, mfc="white", mec="black")
    bx.legend(loc="upper left", fontsize=8, frameon=False)
    bx.set_xticks(range(len(lv)), [f"{f:.1f}" for f in lv])
    bx.xaxis.set_minor_locator(NullLocator())
    bx.set_xlim(-0.5, len(lv) - 0.5)
    style_axes(ax, ylabel=ylabel)
    style_axes(bx, xlabel="Fmax (pN; levels, not to scale)",
               ylabel=f"Shift from {levels[0]:.0f} pN ({unit})")
    ax.tick_params(labelbottom=False)
    for a_ in (ax, bx):
        a_.tick_params(axis="x", labelsize=8)
    # two stacked panels share the x label; the per-axis label check is off
    export_panel(fig, out / "fig_stall_by_fmax", dpi=600, require_axis_labels=False)


def fig_histograms(t, out, xlabel, show_sites, unit):
    """Per force: fraction of that force's tethers in each 30 bin (n falls with
    force, so counts are not comparable between rows; user, 2026-10-08)."""
    n_lv = t.groupby("fmax_pn").tether.nunique()
    levels = [f for f in LEVELS_ALL if n_lv.get(f, 0) >= MIN_TETHERS_PER_LEVEL]
    edges = np.arange(YLIM[0], YLIM[1] + 1, 30)
    hists = {}
    for f in levels:
        x = t[t.fmax_pn == f].stall.to_numpy(float)
        inr = (x >= edges[0]) & (x < edges[-1])
        hists[f] = (x, np.histogram(x[inr], bins=edges)[0] / len(x))
    ymax = max(h.max() for _, h in hists.values())
    fig, axes = plt.subplots(len(levels), 1, figsize=(3.39, 0.8 * len(levels) + 0.7),
                             sharex=True, squeeze=False)
    for ax, f in zip(axes[:, 0], levels):
        x, h = hists[f]
        if show_sites:
            shade(ax, horizontal=False)
        ax.stairs(h, edges, color="black", lw=1.0)
        ax.axvline(np.median(x), ymin=0.035, ymax=0.965, color=PALETTE["dark_grey"], lw=0.9, ls="--")
        raw_trace(ax, np.clip(x, edges[0], edges[-1]), np.full(len(x), -0.1 * ymax), marker="|",
                  filled=True, marker_area=12, color=PALETTE["dark_grey"], edge_width=0.6)
        nover = int((x >= edges[-1]).sum())
        lab = f"{f:.1f} pN\nn = {len(x)}" + (f"\n{nover} > {edges[-1]:.0f}" if nover else "")
        ax.text(1.02, 0.5, lab, transform=ax.transAxes, ha="left", va="center", fontsize=8)
        ax.set_ylim(-0.18 * ymax, 1.12 * ymax)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xlim(edges[0], edges[-1])
    for k, ax in enumerate(axes[:, 0]):
        style_axes(ax, xlabel=xlabel if k == len(levels) - 1 else None)
    fig.supylabel(f"Fraction of tethers per 30-{unit} bin", fontsize=9)
    fig.tight_layout(h_pad=0.15, rect=(0, 0, 0.82, 1))
    export_panel(fig, out / "fig_stall_histograms", dpi=600, require_axis_labels=False)


def fig_601(full, sites, out, xlabel):
    """Full-length tethers: every stall (row = tether), and each tether's
    fraction inside the 601 sites against the fraction the sites cover."""
    reach = 0.9 * NORM_BP
    v = full[(full.stall >= 0) & (full.stall < reach)]
    frac = v.groupby("tether").stall.apply(lambda x: x.apply(inside).mean()).sort_values()
    order = {k: i for i, k in enumerate(frac.index)}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(3.39, 2.6), sharey=True,
                                 gridspec_kw={"width_ratios": [2.4, 1], "wspace": 0.08})
    shade(ax, horizontal=False)
    shown = full[full.tether.isin(order)]          # tethers with a stall in 0..0.9 x 1028
    inr = shown.stall <= YLIM[1]
    raw_trace(ax, shown.stall[inr], shown.tether.map(order)[inr],
              filled=True, marker_area=8, color=PALETTE["dark_grey"])
    if (~inr).any():
        raw_trace(ax, np.full((~inr).sum(), YLIM[1] + 0.025 * (YLIM[1] - YLIM[0])),
                  shown.tether.map(order)[~inr], marker=">", filled=True, marker_area=8,
                  color=PALETTE["dark_grey"], clip_on=False)
    ax.axvline(reach, ymin=0.03, ymax=0.97, color=PALETTE["light_grey"], ls="--", lw=0.8)
    ax.set_xlim(YLIM[0], YLIM[1])
    ax.set_yticks([])
    ax.set_ylim(-1, len(order))
    raw_trace(bx, frac.values, np.arange(len(frac)), filled=True, marker_area=10,
              color=PALETTE["dark_grey"])
    s0 = sites[sites.site_shift_bp == 0].iloc[0]
    bx.axvline(s0.cover, color="black", ls="--", lw=0.8)
    bx.axvspan(s0.ci_lo, s0.ci_hi, color=PALETTE["red"], alpha=0.15, lw=0, zorder=0)
    bx.axvline(s0.fraction_inside, color=PALETTE["red"], lw=0.8)
    bx.set_xlim(-0.05, 1.05)
    style_axes(ax, xlabel=xlabel, ylabel=f"Full-length tethers (n = {len(order)})")
    style_axes(bx, xlabel="Fraction in 601")
    export_panel(fig, out / "fig_601_full_length", dpi=600, require_axis_labels=False)


CAPTIONS = {
    "fig_stall_by_fmax": (
        "Top: stall at the end of each ~10 s plateau for every selected tether at every Fmax "
        "(grey lines join one tether; filled = analysed; open = the tether's 1028 file, excluded "
        "on the 1028 scale; triangles = values above {top}, n given in the legend). Black "
        "diamonds: tether-adjusted level estimates (median-type LAD additive model, anchored on "
        "the median tether intercept). Bottom: shift from {low} pN with 95% percentile bootstrap "
        "CI over tethers (1000 resamples); the open diamond is the reference (0 by definition). "
        "Levels with < 5 tethers are not modelled. Force levels are equally spaced, not to scale."),
    "fig_stall_histograms": (
        "Stall at the end of the plateau, fraction of that force's tethers per 30-{unit} bin "
        "(rows comparable despite falling n). Ticks under each row: individual tethers. Dashed "
        "line: median of all tethers at that force. Values above {top} counted at the right. "
        "Rows need >= 5 tethers."),
    "fig_601_full_length": (
        "Full-length tethers only (span 0.80-1.25 x the dWLC full-unzip span). Left: every stall "
        "(row = tether; triangles = above {top}); grey bands: 601 sites 314-460 / 511-657; grey "
        "dashed: 0.9 x 1028 analysis limit. Right: per-tether fraction of stalls inside the sites "
        "(points), their mean (red line) with 95% bootstrap CI over tethers (red band), and the "
        "fraction of 0-925 the sites cover (black dashed)."),
}


def write_captions(out: Path, unit: str, low: float):
    for name, text in CAPTIONS.items():
        if (out / f"{name}.png").exists():
            (out / f"{name}.txt").write_text(
                text.format(top=f"{YLIM[1]:.0f}", unit=unit, low=f"{low:.0f}") + "\n")


def analyse(st_all: pd.DataFrame, key: str, boot_n: int = N_BOOT) -> dict:
    folder, val, ref, ratio, ylabel = SCALES[key]
    out = OUT / folder
    out.mkdir(parents=True, exist_ok=True)
    st = scale_view(st_all, key)
    t = levels_table(st)
    t = t[t.groupby("tether").fmax_pn.transform("nunique") >= 2]
    refs = ref_levels(st)
    refs = refs[refs.tether.isin(t.tether.unique())]
    t.to_csv(out / "stall_by_tether_level.csv", index=False)
    tm, levels = model_set(t)
    fit_l, base = fit_shift(tm, levels, lad=True)
    fit_s, _ = fit_shift(tm, levels, lad=False)
    b_l = boot_shift(tm, levels, lad=True, n_boot=boot_n)
    rows = []
    for f in sorted(t.fmax_pn.unique()):
        x = t[t.fmax_pn == f].stall
        rows.append({"fmax_pn": f, "n_tethers": len(x), "median": round(x.median(), 0),
                     "q25": round(x.quantile(.25), 0), "q75": round(x.quantile(.75), 0),
                     "in_model": f in levels, "shift_lad": round(fit_l.get(f, np.nan), 0),
                     "lad_lo": round(b_l[f].quantile(.025), 0) if f in b_l else np.nan,
                     "lad_hi": round(b_l[f].quantile(.975), 0) if f in b_l else np.nan,
                     "shift_ls": round(fit_s.get(f, np.nan), 0)})
    shift = pd.DataFrame(rows)
    steps = step_tests(t)
    trend = trend_tests(t)
    cyc = st[st.used & ~st.is_1028_file].pivot_table(index=["tether", "data_id"],
                                                     columns="cycle", values="stall").dropna()
    dcc = cyc[2] - cyc[1]
    within = pd.DataFrame([{"comparison": "C2 - C1, same file and force", "n": len(dcc),
                            "median": round(dcc.median(), 0), "q25": round(dcc.quantile(.25), 0),
                            "q75": round(dcc.quantile(.75), 0), "frac_up": round((dcc > 0).mean(), 2),
                            "wilcoxon_p": wilcoxon(dcc).pvalue}])
    full = t[t.span_ratio.between(*FULL_RATIO)]
    sites = site_tests(full, folder)
    shift.to_csv(out / "stats_levels_and_shift.csv", index=False)
    pd.concat([steps.assign(section="neighbouring forces, same tether"),
               trend.assign(section="per tether"), within.assign(section="within file")],
              ignore_index=True).to_csv(out / "stats_tests.csv", index=False)
    sites.to_csv(out / "stats_601_full_length.csv", index=False)
    unit = "bp" if key == "bp" else "unit"
    fig_trend(t, refs, levels, shift, base, ylabel, out, show_sites=key == "bp",
              unit="bp" if key == "bp" else "1028 scale")
    fig_histograms(t, out, ylabel, show_sites=key == "bp", unit=unit)
    if full.tether.nunique() >= 3 and len(sites):
        fig_601(full, sites, out, ylabel)
    write_captions(out, unit, levels[0])
    print(f"\n===== {folder}: tethers {t.tether.nunique()}, modelled levels {levels}")
    print(shift.to_string(index=False))
    print(steps.to_string(index=False))
    print(trend.to_string(index=False))
    print(within.to_string(index=False))
    print(f"full-length tethers: {full.tether.nunique()}")
    print(sites.to_string(index=False))
    return {"t": t, "shift": shift, "trend": trend, "steps": steps, "sites": sites,
            "within": within, "levels": levels}


def compare_selection(sel_st: pd.DataFrame, unsel_st: pd.DataFrame) -> pd.DataFrame:
    """The same rule-free (dWLC bp) trend in the selected and the unselected
    normalised tethers of the trimmed fields, overall and quality-matched."""
    rows = []
    for name, st in (("selected", sel_st), ("unselected", unsel_st)):
        t = levels_table(scale_view(st, "bp"))
        t = t[t.groupby("tether").fmax_pn.transform("nunique") >= 2]
        span = st.drop_duplicates("tether").set_index("tether").span_A_nm
        for sub, tt in (("all with >= 2 forces", t),
                        ("span >= 150 nm and >= 3 forces",
                         t[t.tether.map(span).ge(150)
                           & t.groupby("tether").fmax_pn.transform("nunique").ge(3)])):
            if tt.tether.nunique() < 5:
                continue
            tr = trend_tests(tt)
            tm, lv = model_set(tt)
            sh = np.nan
            if 18.0 in lv and 39.9 in lv:
                sh = fit_shift(tm, lv, lad=True)[0][39.9]
            rows.append({"set": name, "subset": sub, "n_tethers": tt.tether.nunique(),
                         "frac_highest_above_lowest": tr.frac_up.iloc[0],
                         "median_high_minus_low_bp": tr["median"].iloc[0],
                         "spearman_median": tr["median"].iloc[1],
                         "frac_rho_positive": tr.frac_up.iloc[1],
                         "lad_shift_18_to_39.9_bp": round(sh, 0)})
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    apply_style()
    idx = pd.read_csv(V2 / "norm1028_beads.csv")
    # the user's current selection (66 tethers), column "selected"
    sel = idx[idx.selected.eq("y")]
    fields = set(sel.fov)
    unsel = idx[idx.fov.isin(fields) & ~idx.selected.eq("y") & idx.normalised]
    a = pd.read_csv(DEST / "per_cycle_anchors.csv", low_memory=False)
    a["level_jump_bad"] = a.level_jump_bad.astype(str).eq("True")
    first = a.drop_duplicates(["data_id", "roi_key"])
    print(f"selected tethers: {len(sel)} {sel.fov.value_counts().sort_index().to_dict()}")
    st = collect(sel, first)
    st.to_csv(OUT / "stall_by_cycle.csv", index=False)
    tb = st.drop_duplicates("tether")
    above = st[st.used & (st.data_id < st.ref_file_A)].groupby(["tether", "fmax_pn"]).stall_A.median()
    print(f"plateau cycles {len(st)}; used {int(st.used.sum())}; not used "
          f"{st[~st.used].why_not_used.value_counts().to_dict()}")
    print(f"rule B differs from A in {int((tb.ref_file_A != tb.ref_file_B).sum())} tethers; "
          f"tethers with a lower-force level above 1028 under A: "
          f"{above[above > NORM_BP].reset_index().tether.nunique()}; full-length A "
          f"{int(tb.span_ratio_A.between(*FULL_RATIO).sum())}, B "
          f"{int(tb.span_ratio_B.between(*FULL_RATIO).sum())} of {len(tb)}")
    res = {k: analyse(st, k) for k in SCALES}
    print(f"\nunselected normalised tethers in the trimmed fields: {len(unsel)}")
    ust = collect(unsel, first)
    ust.to_csv(OUT / "stall_by_cycle_unselected.csv", index=False)
    comp = compare_selection(st, ust)
    comp.to_csv(OUT / "selection_comparison.csv", index=False)
    print(comp.to_string(index=False))
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
