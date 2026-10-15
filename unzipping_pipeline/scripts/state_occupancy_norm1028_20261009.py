#!/usr/bin/env python3
"""20261009 dinucleosome: fork-state occupancy vs force and vs time in the plateau (v2, after audit).

User (2026-10-11): "先做状态占据随力和时间变化的图" - the systematic summary of where the
fork is (closed / before N1 / at N1 / at N2 / past N2 / fully open) at every force of the
ladder and at every moment of the 10 s hold.

Data: the user's 96 tethers (selection_96_tethers.csv); every step the bead survived
(ladder_norm1028_20261009.bead_scale: zero = max extension of the 10 pN step, rest-dwell
drift, steps after a loss / fly-away / an unclosed rest dropped). Signal: the display-smoothed
Z (median 5 + mean 5 frames) minus the rest-dwell drift minus the zero, nm. t = 0 is the end of
the 2 mm/s up move; the up move is kept from 10 pN on as negative time (below 10 pN the fork
cannot open and the tracking at 1.6-5 pN is unreliable). data_003 (15 s holds) is cut at 10 s.

SCALE (env): norm1028 (default, the user's: 1028 * (ext - zero) / span) or modelB (absolute
opened bp, 2 x_ss per bp after the handle stretch from 10 pN) as a robustness check - the two
disagree for the 65 / 96 tethers whose 1028 point is < 874 bp in model B. A 0.25 s running
median is applied on the chosen scale.

States (units of the scale; YTrunk_2601_980bp, 601 sites 314-460 / 511-657):
  closed < 100 | before N1 100-280 | at N1 280-520 | at N2 520-690 | past N2 690-950 | full >= 950
v1 split N1 / N2 at 485 (linker middle); the audit found 485 inside a continuous band of >= 1 s
holds (13 / 32 tethers stalling at 445-525 on both sides of it), the first gap is at 520-540
(N2 entry 511), so v2 uses 520. Threshold sensitivity: history_threshold_sensitivity.csv.

Passage (v2): a threshold counts as passed at the first sample at or above it after which the
smoothed position stays >= threshold - 30 for 0.5 s (or to the end of the hold); a single
excursion of a fork parked on a boundary no longer counts.

Left out (v2, audit): (1) every tether's reference step - its 1028 point is the maximum of that
same plateau, so it passes N1 and N2 by construction (61 of 62 plateaus at 39.9 pN); (2) tethers
whose span is < 0.3 x a full 971 bp unzip at their reference force (model B) or whose 10 pN
closed-fork hold moves > 140 units (99th percentile): noise amplified on the 0-1028 scale;
(3) plateaus whose median is > 1150 or < -300 (tracking failure). All are listed in
plateau_states.csv / excluded_tethers.csv.

History per barrier (N1: reached = 280, passed = 520; N2: reached = 520, passed = 690):
  first contact   - the fork never reached the barrier in an earlier step;
  repeat contact  - it reached the barrier in an earlier step and was held there (not passed);
  crossed before  - it passed the barrier in an earlier step.
(v1 pooled the first two as 'first encounter'; 59-75 % of that group at >= 26.9 pN were repeat
contacts, i.e. survivors of lower forces.) Proportions: Wilson 95 % intervals; pooled force
groups: bootstrap over tethers.
Crossing force on the up move: the 10 -> 25 pN part of the ramp is the same at every step, so
crossings are compared as cumulative fractions by 17 / 18 / 20 / 22 / 25 pN on steps >= 25.9 pN
(v1's 'median force at crossing during the move' was truncated at each step's target force).
'No stall' (< 0.25 s from reaching to passing) is given next to the naked-DNA ceiling: the
fraction of 100 -> 280 traversals (no barrier) that take < 0.25 s at the same force.

Outputs (WHZ_20261009/state_occupancy_norm1028/ or state_occupancy_modelB/; v1 in v1_before_audit/):
  fig1_end_state_by_force      state at the end of the hold (last 0.5 s), stacked bars per force
  fig2_occupancy_vs_time       state occupancy vs time, one panel per force
  fig3_crossing                fraction past N1 / past N2 vs time (Kaplan-Meier) and by 10 s vs force
  fig4_history                 P(pass | reached) by history, and the up-move crossing forces
  fig5_barrier_dwell           dwell at N1 / N2 on first contact (Kaplan-Meier)
  *.csv                        end_state_by_force, occupancy_vs_time, crossing_by_force,
                               history_by_force, history_pooled, ramp_crossing_cdf,
                               history_threshold_sensitivity, plateau_states, excluded_tethers

    .venv/bin/python state_occupancy_norm1028_20261009.py [--recollect]
    SCALE=modelB .venv/bin/python state_occupancy_norm1028_20261009.py
    BATCH=20261013 [SCALE=modelB] .venv/bin/python state_occupancy_norm1028_20261009.py   (+Mg2+ batch -> WHZ_20261013/)
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("MODEL", "B")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import process_20261009 as p  # noqa: E402

# BATCH (env): 20261009 (default, no Mg2+) or 20261013 (+Mg2+, same ladder; user's 101 tethers).
# The 1009 pipeline is pointed at the other folder before ladder_norm1028 is imported.
BATCH = os.environ.get("BATCH", "20261009")
if BATCH != "20261009":
    p.RAW = Path(f"/Users/whz/Desktop/DATA MT/unzipping/{BATCH}")
    p.DEST = Path(f"/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces/WHZ_{BATCH}")
import ladder_norm1028_20261009 as L  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

SEL_CSVS = {"20261009": L.SEL_CSV, "20261013": p.DEST / "selection_101_tethers.csv"}

import matplotlib.pyplot as plt  # noqa: E402
from whz_plot import PALETTE, style_axes  # noqa: E402

SCALE = os.environ.get("SCALE", "norm1028")
OUT = p.DEST / ("state_occupancy_norm1028" if SCALE == "norm1028" else "state_occupancy_modelB")
CACHE = p.DEST / "state_occupancy_norm1028" / "_plateaus_v2.pkl"      # nm, shared by both scales
UNIT = "0-1028" if SCALE == "norm1028" else "bp, model B"
NORM = 1028.0
TRUNK_BP = 971.0
EDGES = (100.0, 280.0, 520.0, 690.0, 950.0)
NAMES = ("closed", "before N1", "at N1", "at N2", "past N2", "full")
COLS = (PALETTE["light_grey"], PALETTE["cyan"], PALETTE["blue"], PALETTE["red"], PALETTE["orange"],
        PALETTE["dark_grey"])
NAKED_LO, REACH_N1, PAST_N1, PAST_N2 = EDGES[0], EDGES[1], EDGES[2], EDGES[3]
BARRIERS = {"N1": (REACH_N1, PAST_N1), "N2": (PAST_N1, PAST_N2)}
HISTORIES = ("first contact", "repeat contact", "crossed before")
F_GROUPS = ((18.0, 18.0), (25.9, 27.9), (30.0, 32.2), (34.6, 37.1))
SMOOTH_S = 0.25
SUSTAIN_S = 0.5
HYST = 30.0
HOLD_S = 10.0
END_WIN_S = 0.5
F_RAMP_MIN = 10.0
T_UP_MIN = -1.0
GRID = np.round(np.arange(T_UP_MIN, HOLD_S + 1e-9, 0.05), 3)
N_BOOT = 2000
MIN_N = 5
NOISE_MAX = 140.0          # 99th percentile |x - median| on the 10 pN hold, 0-1028 units
SPAN_RATIO_MIN = 0.3
PLATEAU_OK = (-300.0, 1150.0)
RAMP_F = (16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0)
THR_SENS = (435.0, 485.0, 520.0, 585.0)


# ----------------------------------------------------------------- collect (nm, scale-free)
def collect() -> list[dict]:
    sel = pd.read_csv(SEL_CSVS[BATCH])
    plats = []
    for d, g in sel.groupby("data_id"):
        run = p.load(int(d))
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{int(d):03d}.h5")
        t = run["time_s"]
        for bead in g.bead:
            roi = bead[1:]
            ext, rec = L.bead_scale(run, steps, roi, trk.get(roi, {}))
            zero, span = rec["zero_nm"], rec["span_nm"]
            if not (np.isfinite(zero) and np.isfinite(span) and span > 0):
                continue
            full = float(L.model_ext(np.array([rec["ref_force_pn"]]), TRUNK_BP)[0])
            for st in steps:
                if st["step"] not in rec["usable_steps"] or st["force_pn"] > 40.0:
                    continue
                pl = st["plateau"]
                ix = np.concatenate([st["up"], pl])
                tp = t[ix] - t[pl[0]]
                fx = run["force_pn"][ix]
                keep = ((t[ix] < rec["t_loss"]) & np.isfinite(ext[ix]) & (tp <= HOLD_S)
                        & ((tp >= 0) | (fx >= F_RAMP_MIN)))
                if (keep & (tp >= 0)).sum() < 20:
                    continue
                plats.append({"data_id": int(d), "bead": bead, "step": st["step"],
                              "force_pn": round(st["force_pn"], 1), "t": tp[keep],
                              "e": (ext[ix][keep] - zero).astype(float), "f": fx[keep].astype(float),
                              "hold_s": float(tp[keep][-1]), "span_nm": span, "ref_force_pn": rec["ref_force_pn"],
                              "span_ratio": span / full, "is_ref_step": st["step"] == rec["ref_step"]})
        print(f"data_{int(d):03d}", flush=True)
    return plats


def load(recollect=False) -> list[dict]:
    if CACHE.exists() and not recollect:
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    res = collect()
    with open(CACHE, "wb") as fh:
        pickle.dump(res, fh)
    return res


def prepare(q: dict) -> dict:
    """Position on the chosen scale (raw and 0.25 s running median) and the 0-1028 copy used by
    the quality gate (identical for SCALE=norm1028)."""
    q = dict(q)
    xn = NORM * q["e"] / q["span_nm"]
    xr = xn if SCALE == "norm1028" else L.model_luz(q["e"], q["f"])
    ok = np.isfinite(xr)
    for k in ("t", "e", "f"):
        q[k] = q[k][ok]
    xr, xn = xr[ok], xn[ok]
    fps = (q["t"].size - 1) / (q["t"][-1] - q["t"][0]) if q["t"].size > 1 and q["t"][-1] > q["t"][0] else 40.0
    w = max(3, int(round(SMOOTH_S * fps)))
    sm = lambda v: pd.Series(v).rolling(w, center=True, min_periods=w // 2 + 1).median().bfill().ffill().to_numpy()
    q.update(x_raw=xr, x=sm(xr), xn=sm(xn))
    return q


# ----------------------------------------------------------------- per plateau
def state_of(x, shift=0.0) -> np.ndarray:
    return np.digitize(x, np.asarray(EDGES) + shift)


def on_grid(q) -> np.ndarray:
    j = np.clip(np.searchsorted(q["t"], GRID, side="right") - 1, 0, q["t"].size - 1)
    x = q["x"][j].astype(float)
    x[(GRID > q["hold_s"]) | (GRID < q["t"][0])] = np.nan
    return x


def passage(q, thr) -> tuple[float, bool, float, float]:
    """Sustained first passage: (time, event, force at the sample, force at the sample before)."""
    t, x, f = q["t"], q["x"], q["f"]
    for i in np.flatnonzero(x >= thr):
        w = (t >= t[i]) & (t <= t[i] + SUSTAIN_S)
        if w.sum() >= 3 and np.all(x[w] >= thr - HYST):
            return float(t[i]), True, float(f[i]), float(f[max(i - 1, 0)])
    return q["hold_s"], False, np.nan, np.nan


def end_position(q) -> float:
    return float(np.median(q["x"][q["t"] >= q["hold_s"] - END_WIN_S]))


def km(times: np.ndarray, events: np.ndarray, grid=GRID) -> np.ndarray:
    """Kaplan-Meier 1 - S(t) (checked against brute force in the audit)."""
    ut = np.unique(times[events])
    surv = np.cumprod([1.0 - np.sum((times == u) & events) / np.sum(times >= u) for u in ut])
    j = np.searchsorted(ut, grid, side="right") - 1
    return np.where(j >= 0, 1.0 - surv[np.clip(j, 0, None)] if ut.size else 0.0, 0.0)


def km_median(dw, ev) -> float:
    if dw.size == 0:
        return np.nan
    g = np.round(np.arange(0.0, HOLD_S + 1.5, 0.01), 2)
    j = np.flatnonzero(km(dw, ev, g) >= 0.5)
    return float(g[j[0]]) if j.size else np.inf


def wilson(k: float, n: float, z=1.96) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    ph = k / n
    c = (ph + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)


def tether_boot(df: pd.DataFrame, col: str, rng) -> tuple[float, float]:
    """95 % interval of a pooled proportion, resampling tethers."""
    g = df.groupby(["data_id", "bead"])[col].agg(["sum", "size"]).to_numpy(float)
    if len(g) < 2:
        return np.nan, np.nan
    b = []
    for _ in range(N_BOOT):
        s = g[rng.integers(0, len(g), len(g))]
        b.append(s[:, 0].sum() / s[:, 1].sum())
    return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))


# ----------------------------------------------------------------- gate + per-plateau table
def gate(plats: list[dict]) -> tuple[list[dict], pd.DataFrame]:
    rows = []
    by = {}
    for q in plats:
        by.setdefault((q["data_id"], q["bead"]), []).append(q)
    bad = {}
    for k, qs in by.items():
        q10 = [q for q in qs if q["step"] == 1]
        noise = np.nan
        if q10:
            h = q10[0]["xn"][q10[0]["t"] >= 0]
            noise = float(np.percentile(np.abs(h - np.median(h)), 99))
        why = []
        if qs[0]["span_ratio"] < SPAN_RATIO_MIN:
            why.append(f"span {qs[0]['span_nm']:.0f} nm = {qs[0]['span_ratio']:.2f} x a full unzip")
        if np.isfinite(noise) and noise > NOISE_MAX:
            why.append(f"10 pN hold moves {noise:.0f} units (99th pct)")
        rows.append({"data_id": k[0], "bead": k[1], "span_nm": round(qs[0]["span_nm"], 0),
                     "span_ratio": round(qs[0]["span_ratio"], 2), "noise99_10pN": round(noise, 0),
                     "excluded": bool(why), "why": "; ".join(why)})
        if why:
            bad[k] = "; ".join(why)
    for q in plats:
        med = float(np.median(q["xn"][q["t"] >= 0]))
        q["excluded"] = ("tether: " + bad[(q["data_id"], q["bead"])] if (q["data_id"], q["bead"]) in bad
                         else "reference step (1028 by construction)" if q["is_ref_step"]
                         else f"plateau median {med:.0f} outside {PLATEAU_OK}" if not PLATEAU_OK[0] < med < PLATEAU_OK[1]
                         else "")
    return plats, pd.DataFrame(rows)


def plateau_table(plats: list[dict]) -> pd.DataFrame:
    rows = []
    for q in plats:
        e = end_position(q)
        r = {"data_id": q["data_id"], "bead": q["bead"], "step": q["step"], "force_pn": q["force_pn"],
             "is_ref_step": q["is_ref_step"], "excluded": q["excluded"], "hold_s": round(q["hold_s"], 2),
             "x_first": round(float(np.median(q["x_raw"][:3])), 0), "end_pos": round(e, 0),
             "end_state": NAMES[int(state_of(e))]}
        for nm_, thr in (("naked100", NAKED_LO), ("reach_N1", REACH_N1), ("past_N1", PAST_N1), ("past_N2", PAST_N2),
                         *[(f"past_{int(s)}", s) for s in THR_SENS if s != PAST_N1]):
            tt, ev, fa, fb = passage(q, thr)
            r[f"t_{nm_}"] = round(tt, 3) if ev else np.nan
            r[f"F_{nm_}"] = round(fa, 2) if ev else np.nan
            r[f"Fbefore_{nm_}"] = round(fb, 2) if ev else np.nan
        r["started_past_N1"] = r["x_first"] >= PAST_N1
        r["started_past_N2"] = r["x_first"] >= PAST_N2
        rows.append(r)
    ps = pd.DataFrame(rows).sort_values(["data_id", "bead", "step"]).reset_index(drop=True)
    # history: from every earlier step of the same tether that is not a tracking failure
    valid = ~ps.excluded.str.startswith("plateau median")
    grp = [ps.data_id, ps.bead]
    prior = lambda c: ((ps[c].notna() & valid).groupby(grp).transform(lambda v: v.shift(fill_value=False).cummax())
                       .astype(bool))
    for b, (ent, ex) in BARRIERS.items():
        ce = "t_reach_N1" if b == "N1" else "t_past_N1"
        cx = "t_past_N1" if b == "N1" else "t_past_N2"
        crossed, reached = prior(cx), prior(ce)
        ps[f"{b}_history"] = np.where(crossed, "crossed before", np.where(reached, "repeat contact", "first contact"))
    for s in THR_SENS:
        c = "t_past_N1" if s == PAST_N1 else f"t_past_{int(s)}"
        ps[f"N1_crossed_before_{int(s)}"] = prior(c)
    return ps


# ----------------------------------------------------------------- tables
def tables(plats, ps):
    rng = np.random.default_rng(20261011)
    use = ps[ps.excluded == ""]
    forces = [f for f, g in use.groupby("force_pn") if len(g) >= MIN_N]
    end, occ, cross, hist, pooled, sens = [], [], [], [], [], []
    for f in forces:
        g = use[use.force_pn == f]
        r = {"force_pn": f, "n_tethers": len(g)}
        for nm_ in NAMES:
            k = int((g.end_state == nm_).sum())
            lo, hi = wilson(k, len(g))
            r.update({nm_: round(k / len(g), 3), f"{nm_}_lo": round(lo, 3), f"{nm_}_hi": round(hi, 3)})
        end.append(r)
        qs = [q for q in plats if q["force_pn"] == f and q["excluded"] == ""]
        xg = np.array([on_grid(q) for q in qs])
        st = state_of(xg).astype(float)
        st[np.isnan(xg)] = np.nan
        n_at = np.sum(np.isfinite(st), axis=0)
        for k, nm_ in enumerate(NAMES):
            occ.append(pd.DataFrame({"force_pn": f, "t_s": GRID, "state": nm_, "n_tracked": n_at,
                                     "fraction": np.round(np.sum(st == k, axis=0) / np.maximum(n_at, 1), 4)}))
        c = {"force_pn": f, "n_tethers": len(g)}
        for key, m in (("P_past_N1", g.t_past_N1.notna()), ("P_past_N2", g.t_past_N2.notna())):
            c[key] = round(m.mean(), 3)
            c[f"{key}_lo"], c[f"{key}_hi"] = (round(v, 3) for v in wilson(m.sum(), len(g)))
        g1 = g[g.t_past_N1.notna()]
        k2 = int(g1.t_past_N2.notna().sum())
        c.update(n_past_N1=len(g1), P_past_N2_given_N1=round(k2 / len(g1), 3) if len(g1) else np.nan)
        c["P_past_N2_given_N1_lo"], c["P_past_N2_given_N1_hi"] = (round(v, 3) for v in wilson(k2, len(g1)))
        c["past_N1_during_move"] = round(float((g.t_past_N1 < 0).mean()), 3)
        c["past_N2_during_move"] = round(float((g.t_past_N2 < 0).mean()), 3)
        cross.append(c)
        # naked-DNA ceiling for 'no stall': 100 -> 280 traversal < 0.25 s
        nk = g[g.t_naked100.notna() & g.t_reach_N1.notna() & (g.x_first < NAKED_LO)]
        ceil = float(((nk.t_reach_N1 - nk.t_naked100) < 0.25).mean()) if len(nk) else np.nan
        for b, (ent, ex) in BARRIERS.items():
            ce = "t_reach_N1" if b == "N1" else "t_past_N1"
            cx = "t_past_N1" if b == "N1" else "t_past_N2"
            for hname in HISTORIES:
                h = g[g[ce].notna() & ~g[f"started_past_{b}"] & (g[f"{b}_history"] == hname)]
                if not len(h):
                    continue
                ev = h[cx].notna().to_numpy()
                dw = np.clip(np.where(ev, h[cx] - h[ce], h.hold_s - h[ce]).astype(float), 0, None)
                lo, hi = wilson(ev.sum(), len(h))
                hist.append({"force_pn": f, "barrier": b, "history": hname, "n_reached": len(h),
                             "P_pass_given_reached": round(ev.mean(), 3), "lo": round(lo, 3), "hi": round(hi, 3),
                             "P_pass_during_move": round(float((h[cx] < 0).mean()), 3),
                             "P_no_stall(<0.25s)": round(float(np.mean(ev & (dw < 0.25))), 3),
                             "naked_ceiling_100to280(<0.25s)": round(ceil, 3), "KM_median_dwell_s": km_median(dw, ev)})
        # threshold sensitivity for N1 (pass threshold moved; reach = 280)
        for s in THR_SENS:
            cx = "t_past_N1" if s == PAST_N1 else f"t_past_{int(s)}"
            for before in (False, True):
                h = g[g.t_reach_N1.notna() & (g.x_first < s) & (g[f"N1_crossed_before_{int(s)}"] == before)]
                if len(h):
                    sens.append({"force_pn": f, "N1_pass_threshold": s, "crossed_before": before, "n": len(h),
                                 "P_pass_given_reached": round(float(h[cx].notna().mean()), 3)})
    for lo_f, hi_f in F_GROUPS:
        g = use[(use.force_pn >= lo_f) & (use.force_pn <= hi_f)]
        for b in BARRIERS:
            ce = "t_reach_N1" if b == "N1" else "t_past_N1"
            cx = "t_past_N1" if b == "N1" else "t_past_N2"
            for hname in HISTORIES:
                h = g[g[ce].notna() & ~g[f"started_past_{b}"] & (g[f"{b}_history"] == hname)].copy()
                if not len(h):
                    continue
                h["ok"] = h[cx].notna().astype(float)
                lo, hi = tether_boot(h, "ok", rng)
                pooled.append({"forces": f"{lo_f:g}" if lo_f == hi_f else f"{lo_f:g}-{hi_f:g}", "barrier": b,
                               "history": hname, "n_plateaus": len(h), "n_tethers": h.groupby(["data_id", "bead"]).ngroups,
                               "P_pass_given_reached": round(h.ok.mean(), 3), "tether_boot_lo": round(lo, 3),
                               "tether_boot_hi": round(hi, 3)})
    # cumulative crossing on the common 10 -> 25 pN ramp (steps >= 25.9 pN, crossings before t = 0)
    ramp = []
    r_ = use[use.force_pn >= 25.9]
    sets = {"naked: reach 280 (fork from < 100)": (r_[r_.x_first < NAKED_LO], "reach_N1")}
    for b in BARRIERS:
        cx = "past_N1" if b == "N1" else "past_N2"
        for hname in HISTORIES:
            sets[f"{b} {hname}"] = (r_[(r_[f"{b}_history"] == hname) & ~r_[f"started_past_{b}"]], cx)
    for lab, (h, c) in sets.items():
        if len(h) < MIN_N:
            continue
        fa = h[f"F_{c}"].where(h[f"t_{c}"] < 0)
        ramp.append({"group": lab, "n": len(h), **{f"by_{x:g}pN": round(float((fa <= x).mean()), 3) for x in RAMP_F},
                     "median_dF_sample_gap_pN": round(float((h[f"F_{c}"] - h[f"Fbefore_{c}"]).where(h[f"t_{c}"] < 0).median()), 2)})
    return (forces, pd.DataFrame(end), pd.concat(occ, ignore_index=True), pd.DataFrame(cross), pd.DataFrame(hist),
            pd.DataFrame(pooled), pd.DataFrame(ramp), pd.DataFrame(sens))


# ----------------------------------------------------------------- figures
def _export(fig, name):
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def legend_handles():
    from matplotlib.patches import Patch
    lab = ["closed (< 100)", "before N1 (100-280)", "at N1 (280-520)", "at N2 (520-690)",
           "past N2 (690-950)", "full (>= 950)"]
    return [Patch(facecolor=c, edgecolor="none", label=l_) for c, l_ in zip(COLS, lab)]


def fig_end(end: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(4.6, 2.7))
    x = np.arange(len(end))
    bottom = np.zeros(len(end))
    for nm_, c in zip(NAMES, COLS):
        v = end[nm_].to_numpy()
        ax.bar(x, v, bottom=bottom, color=c, width=0.75, lw=0)
        bottom += v
    for xi, n in zip(x, end.n_tethers):
        ax.text(xi, 1.02, f"{n}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x, [f"{f:g}" for f in end.force_pn], fontsize=7)
    ax.set_ylim(0, 1.09)
    ax.set_xlim(-0.6, len(end) - 0.4)
    style_axes(ax, xlabel="Force (pN)", ylabel="Fraction of tethers at the end of the hold")
    ax.legend(handles=legend_handles()[::-1], fontsize=6, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_title(f"State at 9.5-10 s ({UNIT}; reference steps left out; n above the bars)", fontsize=7, loc="left")
    _export(fig, "fig1_end_state_by_force")


def fig_time(occ: pd.DataFrame, forces):
    fs = [f for f in forces if f >= 18.0]
    nc = 3
    nr = int(np.ceil(len(fs) / nc))
    fig, axes = plt.subplots(nr, nc, figsize=(7.09, 1.85 * nr + 0.5), sharex=True, sharey=True, squeeze=False)
    for ax, f in zip(axes.flat, fs):
        g = occ[occ.force_pn == f]
        base = np.zeros(GRID.size)
        for nm_, c in zip(NAMES, COLS):
            v = g[g.state == nm_].sort_values("t_s").fraction.to_numpy()
            ax.fill_between(GRID, base, base + v, color=c, lw=0, step="post")
            base += v
        ax.axvline(0, color="black", lw=0.6, ls="--")
        ax.set_title(f"{f:g} pN  (n = {int(g.n_tracked.max())})", fontsize=7, loc="left", pad=2)
        ax.set_xlim(T_UP_MIN, HOLD_S)
        ax.set_ylim(0, 1)
        style_axes(ax)
        ax.tick_params(labelsize=6)
    for ax in list(axes.flat)[len(fs):]:
        ax.remove()
    for ax in axes[-1]:
        ax.set_xlabel("Time from reaching the force (s)", fontsize=7)
    for ax in axes[:, 0]:
        ax.set_ylabel("Fraction of tethers", fontsize=7)
    fig.legend(handles=legend_handles(), fontsize=6, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=0.6, w_pad=0.5)
    _export(fig, "fig2_occupancy_vs_time")


def fig_cross(ps: pd.DataFrame, cross: pd.DataFrame, forces):
    use = ps[ps.excluded == ""]
    fs = [f for f in forces if f >= 18.0]
    cmap = plt.get_cmap("viridis")
    col = {f: cmap(i / max(len(fs) - 1, 1)) for i, f in enumerate(fs)}
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.35), gridspec_kw={"width_ratios": [1, 1, 1.05]})
    for ax, (c_, lab) in zip(axes[:2], (("t_past_N1", f"past N1 (>= {PAST_N1:.0f})"), ("t_past_N2", f"past N2 (>= {PAST_N2:.0f})"))):
        for f in fs:
            g = use[(use.force_pn == f) & ~use[f"started_past_{'N1' if 'N1' in c_ else 'N2'}"]]
            ev = g[c_].notna().to_numpy()
            tt = np.where(ev, g[c_].to_numpy(), g.hold_s.to_numpy())
            ax.step(GRID, km(tt, ev), where="post", color=col[f], lw=1.0, label=f"{f:g}")
        ax.axvline(0, color="black", lw=0.6, ls="--")
        ax.set_xlim(T_UP_MIN, HOLD_S)
        ax.set_ylim(0, 1)
        style_axes(ax, xlabel="Time from reaching the force (s)", ylabel=f"Fraction {lab}")
        ax.tick_params(labelsize=6)
    axes[0].legend(title="pN", fontsize=5.5, title_fontsize=6, frameon=False, loc="lower right", ncol=2)
    ax = axes[2]
    xf = cross[cross.force_pn >= 18.0]
    for key, c, m, lab in (("P_past_N1", PALETTE["blue"], "o", "past N1"), ("P_past_N2", PALETTE["red"], "s", "past N2"),
                           ("P_past_N2_given_N1", PALETTE["dark_grey"], "^", "past N2 | past N1")):
        y = xf[key].to_numpy()
        ax.errorbar(xf.force_pn, y, yerr=[y - xf[f"{key}_lo"].to_numpy(), xf[f"{key}_hi"].to_numpy() - y], fmt=m, ms=3.5,
                    mfc="white", mec=c, color=c, lw=0.8, capsize=1.5, label=lab)
    ax.set_ylim(-0.02, 1.02)
    style_axes(ax, xlabel="Force (pN)", ylabel="Fraction by 10 s (Wilson 95 %)")
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    fig.tight_layout(w_pad=1.0)
    _export(fig, "fig3_crossing")


def fig_history(hist: pd.DataFrame, ramp: pd.DataFrame):
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.45), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    hf = hist[hist.force_pn >= 18.0]
    sty = {"first contact": ("o", "white", -0.2), "repeat contact": ("s", "white", 0.0), "crossed before": ("o", None, 0.2)}
    for ax, b, c in ((axes[0], "N1", PALETTE["blue"]), (axes[1], "N2", PALETTE["red"])):
        for hname, (m, mfc, dx) in sty.items():
            h = hf[(hf.barrier == b) & (hf.history == hname) & (hf.n_reached >= MIN_N)]
            y = h.P_pass_given_reached.to_numpy()
            ax.errorbar(h.force_pn + dx, y, yerr=[y - h.lo.to_numpy(), h.hi.to_numpy() - y], fmt=m, ms=3.3,
                        mfc=mfc or c, mec=c, color=c, lw=0.7, capsize=1.3, label=hname)
        ax.set_ylim(-0.02, 1.02)
        style_axes(ax, xlabel="Force (pN)", ylabel=f"P(past {b} by 10 s | reached {b})")
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=5.5, frameon=False, loc="lower right")
        ax.set_title(f"groups with >= {MIN_N} plateaus; Wilson 95 %", fontsize=6, loc="left")
    ax = axes[2]
    lines = {"naked: reach 280 (fork from < 100)": ("black", "-"), "N1 crossed before": (PALETTE["blue"], "-"),
             "N1 first contact": (PALETTE["blue"], ":"), "N1 repeat contact": (PALETTE["blue"], "--"),
             "N2 crossed before": (PALETTE["red"], "-"), "N2 repeat contact": (PALETTE["red"], "--"),
             "N2 first contact": (PALETTE["red"], ":")}
    xs = np.array(RAMP_F)
    for _, r in ramp.iterrows():
        if r.group not in lines:
            continue
        c, ls = lines[r.group]
        ax.plot(xs, [r[f"by_{x:g}pN"] for x in RAMP_F], color=c, ls=ls, lw=1.0, label=f"{r.group} (n={r.n})")
    ax.set_xlim(RAMP_F[0], RAMP_F[-1])
    ax.set_ylim(0, 1)
    style_axes(ax, xlabel="Force on the up move (pN)", ylabel="Fraction passed by that force")
    ax.tick_params(labelsize=6)
    ax.legend(fontsize=4.8, frameon=False, loc="upper left")
    ax.set_title("common 10->25 pN ramp, steps >= 25.9 pN", fontsize=6, loc="left")
    fig.tight_layout(w_pad=0.8)
    _export(fig, "fig4_history")


def fig_dwell(ps: pd.DataFrame):
    use = ps[ps.excluded == ""]
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.4), sharey=True)
    cmap = plt.get_cmap("viridis")
    g_ = np.round(np.arange(0.0, HOLD_S + 1.5, 0.02), 2)
    for ax, b in zip(axes, BARRIERS):
        ce = "t_reach_N1" if b == "N1" else "t_past_N1"
        cx = "t_past_N1" if b == "N1" else "t_past_N2"
        for i, (lo, hi) in enumerate(F_GROUPS):
            h = use[(use.force_pn >= lo) & (use.force_pn <= hi) & use[ce].notna() & ~use[f"started_past_{b}"]
                    & (use[f"{b}_history"] == "first contact")]
            if len(h) < MIN_N:
                continue
            ev = h[cx].notna().to_numpy()
            dw = np.clip(np.where(ev, h[cx] - h[ce], h.hold_s - h[ce]).astype(float), 0, None)
            lab = f"{lo:g}" if lo == hi else f"{lo:g}-{hi:g}"
            ax.step(g_, 1.0 - km(dw, ev, g_), where="post", color=cmap(i / (len(F_GROUPS) - 1)), lw=1.1,
                    label=f"{lab} pN (n = {len(h)})")
        ax.set_xlim(0, HOLD_S + 1)
        ax.set_ylim(0, 1.02)
        style_axes(ax, xlabel=f"Time since the fork reached {b} (s)", ylabel="Fraction still held" if b == "N1" else None)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6, frameon=False, loc="upper right")
        ax.set_title(f"{b}, first contact only", fontsize=7, loc="left")
    fig.tight_layout(w_pad=1.0)
    _export(fig, "fig5_barrier_dwell")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    p.apply_style()
    plats, excl = gate([prepare(q) for q in load("--recollect" in sys.argv)])
    ps = plateau_table(plats)
    forces, end, occ, cross, hist, pooled, ramp, sens = tables(plats, ps)
    ps.to_csv(OUT / "plateau_states.csv", index=False)
    excl.to_csv(OUT / "excluded_tethers.csv", index=False)
    end.to_csv(OUT / "end_state_by_force.csv", index=False)
    occ.to_csv(OUT / "occupancy_vs_time.csv.gz", index=False)
    cross.to_csv(OUT / "crossing_by_force.csv", index=False)
    hist.to_csv(OUT / "history_by_force.csv", index=False)
    pooled.to_csv(OUT / "history_pooled.csv", index=False)
    ramp.to_csv(OUT / "ramp_crossing_cdf.csv", index=False)
    sens.to_csv(OUT / "history_threshold_sensitivity.csv", index=False)
    fig_end(end)
    fig_time(occ, forces)
    fig_cross(ps, cross, forces)
    fig_history(hist, ramp)
    fig_dwell(ps)
    pd.set_option("display.width", 250)
    print(f"scale {SCALE}: {ps.groupby(['data_id', 'bead']).ngroups} tethers, {len(ps)} plateaus; excluded tethers "
          f"{int(excl.excluded.sum())}; plateaus used {int((ps.excluded == '').sum())}; "
          f"left out: {ps.excluded.str.split(':').str[0].replace('', 'used').value_counts().to_dict()}")
    print(end[["force_pn", "n_tethers", *NAMES]].to_string(index=False))
    print(cross.to_string(index=False))
    print(pooled.to_string(index=False))
    print(ramp.to_string(index=False))
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
