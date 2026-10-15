#!/usr/bin/env python3
"""20261009: every CONVINCING small jump of the user's 96 tethers, zoomed in, 0-1028 scale.

User (2026-10-11): "你给我的图太复杂了，先生成1对于所有有确凿小跳的，zoom in！" /
"给出小跳的trace细节". fig1_examples of smalljumps_modelB_20261009.py was 12 hand-picked
panels on the Luz-B scale with 601 bands and dBIC labels; this is the plain version:
one row per plateau - left the whole plateau (zoom window shaded), right the zoom at frame
resolution - and nothing else.

Signal and steps are those of smalljumps_modelB_20261009.py (read only, cached plateaus):
raw Z minus the rest-dwell drift minus the 10 pN zero, single-frame spikes dropped, OP
step fit lambda = 4 from 0.3 s after arrival. Position on the user's scale:
1028 * (ext - zero) / span (span = 1028 point - zero, norm1028_selected_96_peakmap.csv).

Small jump: |step| < 100 on the 0-1028 scale. Convincing ("确凿"), fixed before looking:
  1. both flanking levels >= 10 frames (a level change, not a one- or two-frame excursion);
  2. no camera gap > 0.1 s across the jump (the jump is not hidden between frames);
  3. |step| >= 2 x frame noise sigma (visible in the raw frames, not only in the fit);
  4. abrupt: dBIC >= 6, i.e. two flat levels + one jump beat every continuous
     piecewise-linear curve with <= 3 knots on the two flanking levels (+-3 s) -
     strong evidence against slow creep cut into a step;
  5. not a copy of an overlapping ROI (< 60 px) and not in the same frame as a change
     point of another bead of the file (stage / focus artefact).
Fork jumps: from an open-fork level (50-900 bp B), or from the closed level into the open
range (the fork starts) / from fully open back into it.
False-positive check: the same criteria on the same tethers where no fork can move -
the 10 pN plateaus (closed fork) and >= 18 pN jumps whose BOTH levels are closed (< 50 bp B)
or both fully open (>= 900 bp B) - as a rate per second, scaled to the open-fork time.
(A first version counted 'from < 50 bp' as no fork; its two hits were the fork leaving
the closed level, 27 -> 61 and -19 -> 95 bp B, i.e. real fork events.)

Outputs (WHZ_20261009/smalljumps_modelB/norm1028_scale/zoom/):
  small_jumps_convincing_zoom.pdf (+ page PNGs)   one row per plateau with convincing jump(s)
  small_jumps_convincing.csv                      trace details per jump
  small_jump_candidates_all.csv                   every fork small jump (open fork, or leaving the closed
                                                  level / re-entering from fully open) with pass/fail
  false_positive_check.csv                        criteria applied to the no-fork controls
  small_jumps_audit.csv (written from the audit workflow)   per jump: artefact lens, discrete-signal lens,
                                                  95 % upper bound of a ramp duration, tier, reason
  small_jumps_final_zoom.pdf (+ page PNGs)        after the audit: 'robust' (dark blue) first, then
                                                  'abrupt?' = level change real, abruptness not shown (pale blue)

    .venv/bin/python smalljumps_zoom_norm1028_20261009.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import smalljumps_modelB_20261009 as SJ  # noqa: E402  (read only: cache, OP fit, dBIC)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from whz_plot import raw_trace, fit_line, style_axes  # noqa: E402

D = SJ.p.DEST
IN = SJ.OUT
OUT = IN / "norm1028_scale" / "zoom"
MIN_FRAMES = 10
MAX_GAP_S = 0.1
MIN_SNR = 2.0
MIN_DBIC = 6.0
HALF_WIN_S = 1.5
ROWS_PER_PAGE = 8
C_RAW, C_FIT, C_FWD, C_BACK, C_BIG = "#8C8C8C", "black", "#3C5488", "#00A087", "#E64B35"
C_FWD_PALE = "#8FA3C8"     # level change real, abruptness not shown (audit)
AUDIT_CSV = "small_jumps_audit.csv"


def spans() -> pd.DataFrame:
    s = pd.read_csv(D / "norm1028_selected_96_peakmap.csv")
    s["tether"] = s.data_id.astype(str) + "_" + s.bead
    return s.set_index("tether")[["span_nm", "zero_10pN_max_nm"]]


def refit(q):
    t, y = SJ.search(q)
    sig = SJ.SB.sigma_diff(y)
    return t, y, sig, SJ.fit_levels(t, y, sig, SJ.LAM)


def neighbours() -> dict:
    pr = pd.read_csv(IN / "overlapping_roi_pairs.csv")
    pr = pr[pr.a_present & pr.b_present]
    nb: dict = {}
    for r in pr.itertuples():
        nb.setdefault((r.data_id, r.roi_a), set()).add(r.roi_b)
        nb.setdefault((r.data_id, r.roi_b), set()).add(r.roi_a)
    return nb


def screen(steps: pd.DataFrame, plats, sp: pd.DataFrame) -> pd.DataFrame:
    """Criteria 1-5 on every user-96 change point that is small on the 0-1028 scale."""
    s = steps[(steps.set == "user96") & steps.use & steps.tether.isin(sp.index)].copy()
    s["span_nm"] = s.tether.map(sp.span_nm)
    s = s[s.span_nm > 0]
    s["step_n"] = 1028 * s.step_nm / s.span_nm
    s = s[s.step_n.abs() < 100].copy()
    s["snr"] = s.step_nm.abs() / s.sigma_nm
    s["c1_frames"] = (s.n_before >= MIN_FRAMES) & (s.n_after >= MIN_FRAMES)
    s["c2_gap"] = s.gap_s <= MAX_GAP_S
    s["c3_snr"] = s.snr >= MIN_SNR
    pre = s.c1_frames & s.c2_gap & s.c3_snr
    # 4: dBIC (re-computed for all, also for the steps that are 'big' on the Luz-B scale)
    s["dbic"] = np.nan
    fits = {}
    for i, g in s[pre].groupby("i"):
        q = plats[int(i)]
        t, y, sig, st = fits.setdefault(int(i), refit(q))
        for ix, r in g.iterrows():
            assert abs(q["t_file_s"] + t[int(r.idx)] - r.t_abs_s) < 1e-6, "refit differs from steps_all"
            s.at[ix, "dbic"] = SJ.event_dbic(t, y, st, int(r.k), sig)[0]
    s["c4_abrupt"] = s.dbic >= MIN_DBIC
    # 5: copies of overlapping ROIs and same-frame change points of other beads
    nb = neighbours()
    key = {(q["data_id"], q["bead"], q["step"]): i for i, q in enumerate(plats)}
    allc = steps[["data_id", "bead", "step", "t_abs_s"]]
    grp = {k: g for k, g in allc.groupby(["data_id", "step"])}
    copy, other = [], []
    need = sorted({key[(r.data_id, b, r.step)] for r in s.itertuples()
                   for b in nb.get((r.data_id, r.bead), ()) if (r.data_id, b, r.step) in key})
    ch = SJ.ungated_changes(plats, need)
    for r in s.itertuples():
        hit = False
        for b in nb.get((r.data_id, r.bead), ()):
            j = key.get((r.data_id, b, r.step))
            if j is not None and j in ch:
                hit |= bool(np.any(np.abs(ch[j] - r.t_abs_s) <= r.tol_s))
        copy.append(hit)
        g = grp[(r.data_id, r.step)]
        o = g[g.bead != r.bead].t_abs_s.to_numpy()
        other.append(int(np.sum(np.abs(o - r.t_abs_s) <= r.tol_s)))
    s["copy_of_neighbour"] = copy
    s["n_other_beads_same_frame"] = other
    s["c5_no_artefact"] = ~s.copy_of_neighbour & (s.n_other_beads_same_frame == 0)
    cs = ["c1_frames", "c2_gap", "c3_snr", "c4_abrupt", "c5_no_artefact"]
    s["convincing"] = s[cs].all(axis=1)
    s["fails"] = s[cs].apply(lambda r: " ".join(c.split("_", 1)[1] for c in cs if not r[c]), axis=1)
    lo_, hi_ = SJ.OPEN_LO, SJ.OPEN_HI
    nofork = ((s.from_bp < lo_) & (s.to_bp < lo_)) | ((s.from_bp >= hi_) & (s.to_bp >= hi_))
    s["region"] = np.where(s.force_pn < 18, "10 pN closed fork",
                           np.where(s.open_fork, "open fork",
                                    np.where(nofork, "no fork (closed / fully open)",
                                             "fork leaves closed / re-enters from fully open")))
    return s


def fp_check(s: pd.DataFrame, sp: pd.DataFrame) -> pd.DataFrame:
    pl = pd.read_csv(IN / "plateau_summary.csv")
    pl = pl[pl.set == "user96"].copy()
    pl["tether"] = pl.data_id.astype(str) + "_" + pl.bead
    pl = pl[pl.tether.isin(sp.index)]
    lo = pl.force_pn < 18
    expo = {"open fork": pl.loc[~lo, "exp_open_s"].sum(),
            "no fork (closed / fully open)": (pl.loc[~lo, "exp_closed_s"] + pl.loc[~lo, "exp_full_s"]).sum(),
            "10 pN closed fork": pl.loc[lo, "exp_total_s"].sum()}
    expo["fork leaves closed / re-enters from fully open"] = np.nan   # an event, no exposure of its own
    rows = []
    for reg, e in expo.items():
        g = s[s.region == reg]
        rows.append({"region": reg, "exposure_s": round(e, 1), "n_small_changes": len(g),
                     "n_convincing": int(g.convincing.sum()),
                     "rate_convincing_per_s": round(g.convincing.sum() / e, 5) if e > 0 else np.nan})
    fp = pd.DataFrame(rows)
    ctl = fp[fp.region.isin(["no fork (closed / fully open)", "10 pN closed fork"])]
    n_ctl, e_ctl = int(ctl.n_convincing.sum()), float(ctl.exposure_s.sum())
    rate_ctl = n_ctl / e_ctl
    up95 = (3.0 if n_ctl == 0 else n_ctl + 1.96 * np.sqrt(n_ctl) + 1.0) / e_ctl     # rule of three if 0
    e_open = expo["open fork"]
    n_open = int(fp.loc[fp.region == "open fork", "n_convincing"].iloc[0])
    fp["expected_false_in_open_fork"] = np.nan
    fp.loc[fp.region == "open fork", "expected_false_in_open_fork"] = round(rate_ctl * e_open, 2)
    fp["expected_false_upper95"] = np.nan
    fp.loc[fp.region == "open fork", "expected_false_upper95"] = round(up95 * e_open, 2)
    fp.attrs["summary"] = (f"no-fork controls (both levels closed or both fully open, and the 10 pN plateaus): "
                           f"{n_ctl} convincing in {e_ctl:.0f} s -> {rate_ctl * e_open:.1f} expected false "
                           f"(95 % upper {up95 * e_open:.1f}) in {e_open:.0f} s open fork; observed {n_open}")
    return fp


def big_steps(steps: pd.DataFrame, i: int) -> pd.DataFrame:
    return steps[(steps.i == i) & (steps.cls == "big_fwd")]


def page(groups: list, plats, steps, sp) -> plt.Figure:
    """One row per plateau: left the whole plateau (zoom window shaded), right the zoom
    (+-1.5 s around its convincing jump(s))."""
    fig = plt.figure(figsize=(7.09, 1.18 * ROWS_PER_PAGE + 0.25))
    gs = fig.add_gridspec(ROWS_PER_PAGE, 2, width_ratios=[1.0, 2.1], hspace=0.95, wspace=0.16,
                          left=0.08, right=0.985, top=0.975, bottom=0.05)
    for n, (i, ev) in enumerate(groups):
        r = ev.iloc[0]
        q = plats[int(i)]
        span = sp.span_nm[r.tether]
        dz = q["zero_nm"] - sp.zero_10pN_max_nm[r.tether]          # cache zero vs table zero (< 0.1 nm)
        t, y = q["t"], 1028 * (q["y"] + dz) / span
        ts, ys, sig, st = refit(q)
        lv = 1028 * (st.level_nm.to_numpy() + dz) / span
        lo, hi = max(0.0, ev.t_s.min() - HALF_WIN_S), min(t[-1], ev.t_s.max() + HALF_WIN_S)
        bg = big_steps(steps, i)
        for c, ax in enumerate((fig.add_subplot(gs[n, 0]), fig.add_subplot(gs[n, 1]))):
            w = slice(None) if c == 0 else (t >= lo) & (t <= hi)
            raw_trace(ax, t[w], y[w], color=C_RAW, filled=True, marker_area=1.2 if c == 0 else 3.0,
                      edge_width=0.0, alpha=0.7, rasterized=True)
            xs, yl = [], []
            for k, s_ in enumerate(st.itertuples()):
                a_, b_ = (s_.t0, s_.t1) if c == 0 else (max(s_.t0, lo), min(s_.t1, hi))
                if b_ > a_:
                    xs += [a_, b_, np.nan]
                    yl += [lv[k], lv[k], np.nan]
            fit_line(ax, xs, yl, color=C_FIT, linewidth=0.8)
            for b in bg.itertuples():
                if c == 0 or lo <= b.t_s <= hi:
                    fit_line(ax, [b.t_s, b.t_s], [lv[int(b.k)], lv[int(b.k) + 1]], color=C_BIG, linewidth=1.4)
            for e in ev.itertuples():
                col = C_FWD if e.step_n > 0 else C_BACK
                if getattr(e, "tier", "robust") != "robust":
                    col = C_FWD_PALE
                fit_line(ax, [e.t_s, e.t_s], [lv[int(e.k)], lv[int(e.k) + 1]], color=col, linewidth=2.0 if c else 1.4)
            if c == 0:
                ax.axvspan(lo, hi, color="#EFEFEF", lw=0, zorder=0)
                ax.set_xlim(0, t[-1])
                yy = y[np.isfinite(y)]
                ax.set_ylim(np.percentile(yy, 0.5) - 40, np.percentile(yy, 99.5) + 40)
                ax.set_title(f"D{r.data_id:03d} {r.bead}  {r.force_pn:g} pN", fontsize=6.5, loc="left", pad=2)
            else:
                yw = y[w]
                ax.set_xlim(lo, hi)
                ax.set_ylim(np.percentile(yw, 0.5) - 25, np.percentile(yw, 99.5) + 25)
                txt = []
                for e in ev.itertuples():
                    near = (f"{abs(e.dt_nearest_big_s):.1f} s {'after' if e.dt_nearest_big_s > 0 else 'before'} "
                            f"{'arrival' if e.nearest_big_kind == 'arrival' else 'big jump'}")
                    tag = ""
                    if hasattr(e, "tier"):
                        tag = "robust  " if e.tier == "robust" else "abrupt?  "
                        if np.isfinite(e.tup):
                            tag += f"ramp \u2264{e.tup:.0f} ms  "
                    txt.append(f"{tag}{e.step_n:+.0f} ({e.step_nm:+.0f} nm, {e.snr:.1f}\u03c3)  {e.from_n:.0f}\u2192{e.to_n:.0f}  "
                               f"dwell {e.dwell_before_s:.2f}/{e.dwell_after_s:.2f} s  {near}")
                ax.set_title("\n".join(txt), fontsize=6.0 if len(ev) == 1 else 5.4, loc="left", pad=2)
            last = n == len(groups) - 1
            style_axes(ax, xlabel="Time in plateau (s)" if last else None,
                       ylabel="Position (0-1028)" if c == 0 else None)
            ax.tick_params(labelsize=5.5)
            ax.xaxis.label.set_size(6.5)
            ax.yaxis.label.set_size(6.5)
    return fig


def final_figure(s: pd.DataFrame, plats, steps, sp):
    """After the adversarial audit (small_jumps_audit.csv: artefact lens + discrete-signal lens per
    jump, plus the near-miss review): robust jumps first, then 'level change real, abruptness
    not shown'; rejected jumps are only listed in the csv."""
    f_ = OUT / AUDIT_CSV
    if not f_.exists():
        return
    au = pd.read_csv(f_)
    au = au[au.tier != "rejected"]
    ev = s.merge(au[["i", "k", "tier", "transition_upper95_ms", "source"]], on=["i", "k"])
    ev["tup"] = ev.transition_upper95_ms.astype(float)
    ev["rank"] = (ev.tier != "robust").astype(int)
    order = ev.groupby("i")["rank"].min().sort_values(kind="stable")
    groups = [(i, ev[ev.i == i].sort_values("t_s")) for i in order.index]
    for g_ in OUT.glob("small_jumps_final_zoom_p*.png"):
        g_.unlink()
    with PdfPages(OUT / "small_jumps_final_zoom.pdf") as pdf:
        for pg in range(0, len(groups), ROWS_PER_PAGE):
            fig = page(groups[pg:pg + ROWS_PER_PAGE], plats, steps, sp)
            pdf.savefig(fig, dpi=300)
            fig.savefig(OUT / f"small_jumps_final_zoom_p{pg // ROWS_PER_PAGE + 1}.png", dpi=200)
            plt.close(fig)
    print(f"final figure: {int((ev.tier == 'robust').sum())} robust + {int((ev.tier != 'robust').sum())} "
          f"'abrupt?' jumps on {len(groups)} plateaus")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    SJ.p.apply_style()
    data = SJ.load_cache()
    plats = data["plateaus"]
    sp = spans()
    steps = pd.read_csv(IN / "steps_all_flagged.csv.gz")
    s = screen(steps, plats, sp)
    s["from_n"] = 1028 * s.from_nm / s.span_nm      # same zero as the table (checked below)
    s["to_n"] = 1028 * s.to_nm / s.span_nm
    dz = np.array([plats[int(i)]["zero_nm"] for i in s.i]) - s.tether.map(sp.zero_10pN_max_nm).to_numpy()
    assert np.nanmax(np.abs(dz)) < 0.1, f"cached zero differs from the table by up to {np.nanmax(np.abs(dz)):.1f} nm"
    fp = fp_check(s, sp)
    cand = s[s.region.isin(["open fork", "fork leaves closed / re-enters from fully open"])].sort_values(
        ["data_id", "bead", "force_pn", "t_s"])
    good = cand[cand.convincing]
    s["t_in_file_s"] = s.t_abs_s          # from the file's first frame (H5 'Time (s)' minus its first value)
    cand = s.loc[cand.index]
    good = cand[cand.convincing]
    cols = ["i", "k", "data_id", "bead", "step", "force_pn", "t_s", "t_in_file_s", "step_n", "step_nm", "from_n", "to_n", "snr", "sigma_nm",
            "fps_median", "dwell_before_s", "dwell_after_s", "n_before", "n_after", "before_censored",
            "after_censored", "gap_s", "dbic", "region", "dt_nearest_big_s", "nearest_big_kind", "wclass",
            "has_overlap_neighbour", "copy_of_neighbour", "n_other_beads_same_frame", "convincing", "fails"]
    good[cols].round(3).to_csv(OUT / "small_jumps_convincing.csv", index=False)
    cand[cols].round(3).to_csv(OUT / "small_jump_candidates_all.csv", index=False)
    fp.to_csv(OUT / "false_positive_check.csv", index=False)
    groups = [(i, g.sort_values("t_s")) for i, g in good.groupby("i", sort=False)]
    for f_ in OUT.glob("small_jumps_convincing_zoom_p*.png"):
        f_.unlink()
    with PdfPages(OUT / "small_jumps_convincing_zoom.pdf") as pdf:
        for pg in range(0, len(groups), ROWS_PER_PAGE):
            fig = page(groups[pg:pg + ROWS_PER_PAGE], plats, steps, sp)
            pdf.savefig(fig, dpi=300)
            fig.savefig(OUT / f"small_jumps_convincing_zoom_p{pg // ROWS_PER_PAGE + 1}.png", dpi=200)
            plt.close(fig)
    (OUT / "README.txt").write_text(
        "Every convincing small jump (|step| < 100 on your 0-1028 scale) of your 96 tethers, >= 18 pN, "
        "fork jumps (from an open-fork level, or the fork leaving the closed level). One row per plateau. Left: the whole plateau, grey zoom window. Right: +-1.5 s zoom. Grey dots = "
        "raw frames (spikes > 6 sigma dropped), black = step fit (optimal partitioning, lambda 4, from 0.3 s), "
        "blue / green = the small forward / backward jump, red = big jumps (>= 100 bp B). Title right: size on "
        "0-1028 (nm, multiples of the frame noise), level before -> after, dwell of the two levels, time to the "
        "nearest big jump or the arrival.\nCriteria (fixed before looking): both levels >= 10 frames; no camera gap "
        "> 0.1 s at the jump; |step| >= 2 x frame noise; abrupt (dBIC >= 6 vs any continuous curve with <= 3 "
        "knots on the two levels, +-3 s); no copy of an overlapping ROI, no change point of another bead in the "
        "same frame.\nFalse-positive check (false_positive_check.csv): " + fp.attrs["summary"] + "\n")
    final_figure(s, plats, steps, sp)
    print(fp.to_string(index=False))
    print(fp.attrs["summary"])
    print(cand.groupby("fails").size().sort_values(ascending=False).head(15).to_string())
    print(good[["data_id", "bead", "force_pn", "t_s", "step_n", "step_nm", "snr", "dwell_before_s",
                "dwell_after_s", "dbic", "dt_nearest_big_s", "wclass"]].round(2).to_string(index=False))
    print(f"{len(good)} convincing on {good.data_id.astype(str).add(good.bead).nunique()} tethers -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
