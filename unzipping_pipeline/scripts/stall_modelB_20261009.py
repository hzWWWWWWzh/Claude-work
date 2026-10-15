#!/usr/bin/env python3
"""20261009 ladder: stall statistics in Luz model B, on the user's 96 selected tethers.

User, 2026-10-11: "用模型B重新做停顿统计".

Positions: Luz_B = (extension - zero - 1800 bp handle stretch from 10 pN) / (2 x_ss(F)),
absolute bp from the 10 pN zero (ladder_norm1028_20261009: zero = max extension of
the 10 pN step, drift through the 1.6 pN rests). Model A (2 x_ss - x_ds) is kept
as a column for comparison. Only the steps the bead survived are used (the
ladder_norm1028 rule: a step after a loss, a fly-away or a rest that did not come
back closed is dropped, and so is every later step).

Stall per step = median Luz_B over 9-10 s into the plateau (the last second of a
10 s hold; data_003's 15 s holds are read at the same time), used if the hold
lasted >= 9 s, the bead is tracked through the window and the fork is stationary
there (extension IQR <= 75 nm, ~68 bp B at 30 pN). The 10 pN step is the zero
reference and is not a level in the statistics.

Statistics (as on 20261004, functions from stall_selected_1028_20261004.py): per
force the median and the fractions past each 601 boundary; a median-type (LAD)
additive model stall = tether + force level with 95 % CIs from resampling tethers;
neighbouring steps of the ladder (Wilcoxon, Bonferroni); per tether highest vs
lowest force and Spearman rho; stalls inside the 601 sites against the fraction of
the range they cover, tether as unit, with both sites shifted together as a control.
Selection caveat: the same trend tests on the unselected clean tethers of the same
17 files (tracked or lost, no LUT level-jump artefact), and on a quality-matched
subset (span >= 150 nm, >= 3 forces).

Outputs: WHZ_20261009/stall_modelB/ (each figure with a caption .txt, README.md).

    .venv/bin/python stall_modelB_20261009.py
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
import stall_selected_1028_20261004 as S  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402

OUT = p.DEST / "stall_modelB"
STALL_T0, STALL_T1 = 9.0, 10.0
MAX_IQR_NM = 75.0
SITES = {"N1 entry": 314, "N1 exit": 460, "N2 entry": 511, "N2 exit": 657}
FULL_OPEN_BP = 0.92 * 971


def luz(e_rel, f, model):
    from reprocess_26_core_sop import HANDLE_BP, x_ds, x_ss
    g = 2 * x_ss(f) if model == "B" else 2 * x_ss(f) - x_ds(f)
    return (e_rel - HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, 10.0)))) / g


def collect() -> pd.DataFrame:
    sel = pd.read_csv(L.SEL_CSV)
    chosen = set(zip(sel.data_id.astype(int), sel.bead))
    rows = []
    for d in sorted(sel.data_id.unique()):
        run = p.load(int(d))
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{int(d):03d}.h5")
        t = run["time_s"]
        for roi in run["rois"]:
            tq = trk.get(roi, {})
            bead = f"R{roi}"
            picked = (int(d), bead) in chosen
            if not picked and tq.get("status") not in ("tracked", "lost"):
                continue                        # unselected comparison: clean tethers only
            ext, rec = L.bead_scale(run, steps, roi, tq)
            if not np.isfinite(rec["zero_nm"]):
                continue
            span = rec["span_nm"]
            if np.isfinite(rec["end_nm"]):            # the 1028 point in Luz_B (bp)
                rec["end_luz"] = float(luz(np.array([rec["end_nm"] - rec["zero_nm"]]),
                                           np.array([rec["ref_force_pn"]]), "B")[0])
            for st in steps:
                if st["step"] == 1 or st["step"] not in rec["usable_steps"]:
                    continue
                pl = st["plateau"]
                tp = t[pl] - t[pl[0]]
                hold = float(tp[-1])
                win = (tp >= STALL_T0) & (tp <= STALL_T1) & np.isfinite(ext[pl]) & (t[pl] < rec["t_loss"])
                reached = hold >= STALL_T1 - 0.1 and win.sum() >= 3
                e = ext[pl][win] - rec["zero_nm"]
                f = np.full(e.size, st["force_pn"])
                iqr = float(np.subtract(*np.percentile(e, [75, 25]))) if e.size >= 3 else np.nan
                why = ("hold < 9 s" if hold < 9.0 else
                       "not tracked through 9-10 s" if not reached else
                       "not stationary (IQR > 75 nm)" if not iqr <= MAX_IQR_NM else "")
                lb = float(np.median(luz(e, f, "B"))) if e.size else np.nan
                la = float(np.median(luz(e, f, "A"))) if e.size else np.nan
                rows.append({"data_id": int(d), "bead": bead, "tether": f"D{int(d):03d}_{bead}",
                             "fov": f"D{int(d):03d}", "selected": picked, "tracking": tq.get("status", "?"),
                             "step": st["step"], "fmax_pn": round(st["force_pn"], 1), "hold_s": round(hold, 2),
                             "stall_nm": round(float(np.median(e)), 1) if e.size else np.nan,
                             "stall_B": round(lb, 1), "stall_A": round(la, 1),
                             "stall_1028": round(1028 * float(np.median(e)) / span, 1)
                             if e.size and np.isfinite(span) and span > 0 else np.nan,
                             "span_nm": round(span, 1), "end_luzB_bp": round(rec.get("end_luz", np.nan), 1)
                             if "end_luz" in rec else np.nan,
                             "is_ref_step": st["step"] == rec["ref_step"], "end_iqr_nm": round(iqr, 1),
                             "used": why == "", "why_not_used": why})
        print(f"data_{int(d):03d} done", flush=True)
    return pd.DataFrame(rows)


def per_force(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for f, g in t.groupby("fmax_pn"):
        r = {"fmax_pn": f, "n_tethers": g.tether.nunique(), "median_bp": round(g.stall.median(), 0),
             "q25": round(g.stall.quantile(.25), 0), "q75": round(g.stall.quantile(.75), 0)}
        for k, v in SITES.items():
            r[f"past {k} ({v})"] = round((g.stall > v).mean(), 2)
        r[f"fully open (>= {FULL_OPEN_BP:.0f})"] = round((g.stall >= FULL_OPEN_BP).mean(), 2)
        rows.append(r)
    return pd.DataFrame(rows)


def fig_past_sites(pf: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    cols = [c for c in pf.columns if c.startswith("past") or c.startswith("fully")]
    pal = [PALETTE["blue"], PALETTE["cyan"], PALETTE["orange"], PALETTE["red"], PALETTE["black"]]
    for c, col in zip(cols, pal):
        ax.plot(pf.fmax_pn, pf[c], "o-", ms=3.5, lw=0.9, color=col, label=c)
    ax.set_ylim(0, 1)
    style_axes(ax, xlabel="Force (pN)", ylabel="Fraction of tethers")
    ax.legend(fontsize=5.5, frameon=False, loc="upper left")
    for ext_ in ("pdf", "png"):
        fig.savefig(out / f"fig_fraction_past_sites.{ext_}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    (out / "fig_fraction_past_sites.txt").write_text(
        "Fraction of the selected tethers whose stall (Luz model B, 9-10 s into the hold) lies past each "
        "601 boundary, per force. Lines join the per-force fractions; each point is one force level.\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p.apply_style()
    st = collect()
    st.to_csv(OUT / "stall_by_step.csv", index=False)
    S.OUT = OUT                                      # the 20261004 figure helpers write here
    S.LEVELS_ALL = [18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9, 46.0]
    S.CAPTIONS = {
        "fig_stall_by_fmax": (
            "20261009 ladder, your selected tethers. Top: stall (median Luz, 9-10 s into the 10 s hold) "
            "of every tether at every force step it survived (grey lines join one tether; triangles = above "
            "{top}). Black diamonds: tether-adjusted level estimates (median-type LAD additive model, "
            "anchored on the median tether intercept). Bottom: shift from {low} pN with 95% percentile "
            "bootstrap CI over tethers (1000 resamples). Levels with < 5 tethers are not modelled; force "
            "levels are equally spaced, not to scale. Grey bands: 601 sites 314-460 / 511-657 bp."),
        "fig_stall_histograms": (
            "Stall per force step, fraction of that force's tethers per 30-{unit} bin (n falls with force). "
            "Ticks: individual tethers; dashed: median. Grey bands: 601 sites. Rows need >= 5 tethers."),
    }
    res = {}
    for model, col, sub in (("B", "stall_B", "luz_modelB"), ("A", "stall_A", "luz_modelA_for_comparison"),
                            ("1028", "stall_1028", "norm1028_scale")):
        out = OUT / sub
        out.mkdir(exist_ok=True)
        u = st[st.selected & st.used].assign(stall=lambda d: d[col], cycle=1, span_ratio=np.nan)
        if model == "1028":            # its own reference step reads ~1028 by construction
            u = u[~u.is_ref_step & u.stall.notna()]
        t = (u.groupby(["tether", "fov", "bead", "fmax_pn"], as_index=False)
              .agg(stall=("stall", "median"), n_cycles=("cycle", "size"), span_ratio=("span_ratio", "first")))
        t = t[t.groupby("tether").fmax_pn.transform("nunique") >= 2]
        tm, levels = S.model_set(t)
        fit_l, base = S.fit_shift(tm, levels, lad=True)
        fit_s, _ = S.fit_shift(tm, levels, lad=False)
        b_l = S.boot_shift(tm, levels, lad=True)
        pf = per_force(t)
        pf["in_model"] = pf.fmax_pn.isin(levels)
        pf["shift_lad"] = pf.fmax_pn.map(lambda f: round(fit_l.get(f, np.nan), 0))
        pf["lad_lo"] = pf.fmax_pn.map(lambda f: round(b_l[f].quantile(.025), 0) if f in b_l else np.nan)
        pf["lad_hi"] = pf.fmax_pn.map(lambda f: round(b_l[f].quantile(.975), 0) if f in b_l else np.nan)
        pf["shift_ls"] = pf.fmax_pn.map(lambda f: round(fit_s.get(f, np.nan), 0))
        steps = S.step_tests(t)
        trend = S.trend_tests(t)
        sites = S.site_tests(t, f"luz_model{model}")
        full = t[t.tether.isin(st[st.end_luzB_bp >= 0.9 * 971].tether)] if model != "1028" else t.iloc[:0]
        sites_full = S.site_tests(full, f"luz_model{model}, full-open tethers") if full.tether.nunique() >= 3 else pd.DataFrame()
        t.to_csv(out / "stall_by_tether_level.csv", index=False)
        pf.to_csv(out / "stats_levels_and_shift.csv", index=False)
        pd.concat([steps.assign(section="neighbouring forces, same tether"),
                   trend.assign(section="per tether")], ignore_index=True).to_csv(out / "stats_tests.csv", index=False)
        pd.concat([sites.assign(subset="all selected tethers"),
                   sites_full.assign(subset="tethers that open fully (end >= 0.9 x 971 bp B)")],
                  ignore_index=True).to_csv(out / "stats_601.csv", index=False)
        ylab = "Stall, 0-1028 scale" if model == "1028" else f"Stall, Luz model {model} (bp)"
        S.fig_trend(t, pd.DataFrame(columns=["tether", "fmax_pn", "stall"]), levels, pf.rename(
            columns={"median_bp": "median"}), base, ylab, out, show_sites=True, unit="bp")
        S.fig_histograms(t, out, ylab, show_sites=True, unit="bp")
        S.write_captions(out, "bp", levels[0])
        if model in ("B", "1028"):
            fig_past_sites(pf, out)
        res[model] = (t, pf, steps, trend, sites, sites_full, levels)
        print(f"\n===== model {model}: {t.tether.nunique()} tethers, levels {levels}")
        print(pf.to_string(index=False))
        print(steps.to_string(index=False))
        print(trend.to_string(index=False))
        print(sites.to_string(index=False))
        if len(sites_full):
            print(sites_full.to_string(index=False))
    # selection caveat (model B)
    rows = []
    for name, mask in (("selected", st.selected), ("unselected clean", ~st.selected)):
        u = st[mask & st.used].assign(stall=lambda d: d.stall_B)
        t = u.groupby(["tether", "fov", "fmax_pn"], as_index=False).agg(stall=("stall", "median"))
        t = t[t.groupby("tether").fmax_pn.transform("nunique") >= 2]
        span = st.drop_duplicates("tether").set_index("tether").span_nm
        for sub, tt in (("all with >= 2 forces", t),
                        ("span >= 150 nm and >= 3 forces",
                         t[t.tether.map(span).ge(150) & t.groupby("tether").fmax_pn.transform("nunique").ge(3)])):
            if tt.tether.nunique() < 5:
                continue
            tr = S.trend_tests(tt)
            tm, lv = S.model_set(tt)
            sh = S.fit_shift(tm, lv, lad=True)[0].get(39.9, np.nan) if 18.0 in lv and 39.9 in lv else np.nan
            rows.append({"set": name, "subset": sub, "n_tethers": tt.tether.nunique(),
                         "frac_highest_above_lowest": tr.frac_up.iloc[0],
                         "median_high_minus_low_bp": tr["median"].iloc[0],
                         "spearman_median": tr["median"].iloc[1], "frac_rho_positive": tr.frac_up.iloc[1],
                         "lad_shift_18_to_39.9_bp": round(sh, 0)})
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "selection_comparison.csv", index=False)
    print(comp.to_string(index=False))
    print(f"steps: {len(st)}; selected used {int((st.selected & st.used).sum())}; "
          f"why not used: {st[st.selected & ~st.used].why_not_used.value_counts().to_dict()}")
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
