#!/usr/bin/env python3
"""20261009: brief jumps between two long stalls (user, 2026-10-11:
"注意在两个主要暂停（长时间stalling）中间可能会有一些短暂的跳变，我认为他们很关键。
请你注意他们，着重分析，他们可能是核小体内圈外圈。").

Batch 20261009 ONLY (YTrunk_2601_980bp, 971 bp trunk, N1 601 314-460 dyad 387,
N2 601 511-657 dyad 584). Every titration file holds the whole step ladder;
each 10 s plateau starts from a fork re-zipped during the 15 s rest at
1.6 pN, so every plateau is an independent pull.

Upstream (read-only imports; nothing upstream is edited or overwritten):
  process_20261009      load, segment_steps, titration_files
  fec_all_20261009      smoothed_nm, rest_drift (line through the last 10 s of
                        each 1.6 pN rest), max_extension (0 point = maximum
                        extension of the 10 pN step of cycle 1)
  dwell_steps_20261009  segment / merge: the BIC binary segmentation (BS) that
                        is validated here against the tuned detector
  tracking_quality      per-ROI status / t_loss / level jump

Signal. Each plateau at the FULL frame rate: raw Z minus the rest-dwell drift
line minus the 10 pN zero (nm), no smoothing. Single-frame tracking spikes
(> 6 robust sigma off a +-3 frame running median) are dropped and counted.

Per-bp models (both reported):
  A (SOP dWLC = fec_all_20261009.model_luz): gain 2 x_ss - x_ds per opened bp;
  B (closed trunk not in the force path):    gain 2 x_ss per opened bp
    (B = A / 1.44-1.49 at 18-40 pN).
Absolute Luz = (extension - 1800 bp handle stretch from 10 pN) / gain.

Tethers / plateaus.
  (a) 'excellent': the 39 tethers of select_tethers_20261009 (selected);
  (b) 'all clean': every tracked tether (lost: only complete plateaus before
      the loss) with a finite 10 pN zero.
  Both sets pass the same per-plateau gate: frame noise 5-60 nm (diff MAD),
  LUT flicker <= 1 % frames > 4 sigma off a 9-frame running median, <= 2 %
  identical values (clipped tracking), and every level of the reference fit
  inside -150..1150 bp (model A; a level outside = another LUT level).
  10 pN plateaus (fork closed) are kept as the EXPERIMENTAL null.

Detector. Exact optimal partitioning (OP, PELT pruned) of the extension (nm)
minimising SSE/sigma^2 + lambda ln(n) per change point, levels >= 3 frames,
sigma from the MAD of first differences. BS (dwell_steps_20261009 as is:
PENALTY 3, MIN_FRAMES 4, MERGE 12 bp) is run on the same traces. lambda is
the smallest of 2..6 whose NOISE nulls (N1, N2) give <= 10 % of the real
count of brief events between long stalls (-> lambda = 4).

Validation (matched to every gated plateau: own time grid, own noise):
  N1 white   long-stall skeleton + Gaussian AR(1) noise (own sigma, own lag-1);
  N2 boot    long-stall skeleton + circular block bootstrap (10 frames) of
             its own residuals (keeps ACF to 10 lags and the real tails);
  N3 creep   each long stall its own line + ramps LINEAR IN TIME between long
             stalls over the interval the brief levels occupied (no discrete
             brief level) + bootstrap residuals; answers 'discrete level or
             finite-speed move' through BIC(ramp) - BIC(steps) (stage disc);
  10 pN      real closed-fork plateaus (tracking artefacts, no fork events).
  Injection: forward intermediates (0 -> A -> 2A) and reversible excursions
  (0 -> A -> 0) of 0.02-0.9 s, A = 20-200 nm, sigma 25 / 36 / 50 nm, on REAL
  plateau time grids grouped by median frame rate (100 / 50 / <=33 fps). The
  camera clock is irregular (frame spacing 0.01-0.19 s inside one plateau),
  so a 0.05 s event falls into a frame gap 8 / 21 / 36 % of the time.

Events. Long stall = a level held >= LONG_S = 1 s (the dwell distribution is
continuous, so 1 s is a convention: a 1 s level is located to 3-7 nm; the
counts are re-run at 0.5 and 2 s). Two long stalls are different main stalls
only if they differ by > MOVE_MIN_NM = 40 nm (10 pN control). A brief level
(< 1 s) with a long stall before AND after it on the same plateau is:
  intermediate_fwd  between the two long levels, the second one higher;
  intermediate_back between them, the second one lower;
  excursion_fwd / _back above / below both long levels; 'returns' when the
                    two long levels agree within max(40 nm, 3 SE).
Brief levels before the first / after the last long stall are 'outside'.
Every transition with intermediates gets BIC(ramp) - BIC(steps) (> 6: the
discrete level(s) beat a continuous ramp and a single jump). Every brief
event gets a tracking check (X, Y, LUT amplitude A co-change > 6 robust sigma;
>= 150 nm returning in <= 5 frames) and 'at_closed_fork' (both main stalls
< 100 bp, model A: N1 not reached).

Geometry: detector-aware uniform null (each real transition re-simulated on
its own clock with its intermediate at a uniform position, own residual noise
and durations, then detected); u-shuffle for absolute positions; comb
statistics C_P = <cos(2 pi d1 / P)> for P = 40 / 80 bp with the power of a
'ladder' H1 (intermediate at k x P bp after the preceding stall, re-detected),
also with a 10 % per-tether scale error; excursion amplitudes vs N2 / N3
false events.

Outputs (Proecessed Traces/WHZ_20261009/stall_structure_s549/briefjumps/):
  plateaus_qc.csv                       gate per plateau (noise, ACF, flicker, LUT)
  detector_null_scan_*.csv(.gz)         FP counts per lambda / null; BS for comparison
  detector_injection_recall.csv, detector_detection_limits.csv
  events_between_long.csv               every brief event between long stalls
  events_brief_all.csv, long_stalls.csv, transitions_long_to_long.csv, plateau_counts.csv
  transitions_with_intermediates_ladder.csv   level ladders of forward moves
  discreteness_*.csv(.gz), geometry_*.csv(.gz), rate_by_force.csv,
  long_threshold_sensitivity.csv, stats.json
  fig1_detector_validation  fig2_dwell_threshold  fig3_examples (+ _events.csv)
  fig4_relative_positions   fig5_absolute_positions  fig6_rate_vs_force

    .venv/bin/python stall_structure_briefjumps_20261009.py [--recollect] [--stage all|prep|scan|inject|events|geom|figs]
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pickle
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter

import process_20261009 as p  # noqa: E402  (sets Agg)
import matplotlib.pyplot as plt  # noqa: E402
import fec_all_20261009 as F  # noqa: E402
import dwell_steps_20261009 as DS  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402
from reprocess_26_core_sop import HANDLE_BP, x_ds, x_ss  # noqa: E402

KEY = "briefjumps"

# TETHER_SET=<csv with data_id, bead>: use that list (the user's hand selection,
# selection_96_tethers.csv) in place of the 39 'excellent' tethers; output then
# goes to stall_structure_s549_user96/ ("excellent" in tables/figures = that list)
import os as _os
TETHER_SET = _os.environ.get("TETHER_SET")


def _user_set():
    import pandas as _pd
    s = _pd.read_csv(TETHER_SET)
    return set(zip(s.data_id.astype(int), s.bead.astype(str)))

OUT = p.DEST / ("stall_structure_s549_user96" if TETHER_SET else "stall_structure_s549") / KEY
CACHE = OUT / "_cache_plateaus.pkl"
SPIKE_NSIG, SPIKE_HALF = 6.0, 3
ZERO_F = F.ZERO_F_PN

# detector
LAMBDAS = (2.0, 3.0, 4.0, 5.0, 6.0)
LAM_REF = 3.0                 # reference fit for gating, skeletons and residuals
MMIN = 3
LONG_S = 1.0
LONG_SENS = (0.5, 2.0)
# Two long stalls are the SAME level unless they differ by > MOVE_MIN_NM (or 3 SE):
# on the closed-fork 10 pN plateaus long-to-long level shifts of 15-40 nm occur
# at 0.11 per plateau (drift / tracking, no fork motion) and none above 60 nm,
# so a 'move' between main stalls must exceed 40 nm (~53 bp A / ~36 bp B at 30 pN).
MOVE_MIN_NM = 40.0
DBIC_STRONG = 6.0
# plateau gate
SIG_MIN_NM, SIG_MAX_NM = 5.0, 60.0
FLICKER_MAX = 0.01
REPEAT_MAX = 0.02
GUARD_BP = (-150.0, 1150.0)
# nulls
BLOCK = 10
N_NULL = 3
N_GEOM = 30
FP_TARGET_FRAC = 0.10         # tuned lambda: matched-null false events <= 10 % of the real count

N1, N2 = (314.0, 460.0), (511.0, 657.0)
DYADS = (387.0, 584.0)
LANDMARKS = {"N1 entry": 314, "N1 dyad-40": 347, "N1 dyad": 387, "N1 dyad+40": 427, "N1 exit": 460,
             "N2 entry": 511, "N2 dyad-40": 544, "N2 dyad": 584, "N2 dyad+40": 624, "N2 exit": 657}
LM = np.array(sorted(LANDMARKS.values()), float)
CATS = ("intermediate_fwd", "intermediate_back", "excursion_fwd", "excursion_back")


# ----------------------------------------------------------------- models
def g_a(f):
    f = np.atleast_1d(np.asarray(f, float))
    return 2.0 * x_ss(f) - x_ds(f)


def g_b(f):
    f = np.atleast_1d(np.asarray(f, float))
    return 2.0 * x_ss(f)


def gains(f) -> tuple[float, float]:
    return float(g_a(f)[0]), float(g_b(f)[0])


def handle_nm(f):
    f = np.atleast_1d(np.asarray(f, float))
    return HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, ZERO_F)))


def abs_bp(nm, f, model):
    g = g_a(f) if model == "A" else g_b(f)
    return (np.asarray(nm, float) - handle_nm(f)) / g


def turn_region(x: float) -> str:
    """Position relative to the 601 map (bp): outer entry, inner, outer exit, linker, flank."""
    for name, (lo, hi), dy in (("N1", N1, DYADS[0]), ("N2", N2, DYADS[1])):
        if lo <= x < dy - 40:
            return f"{name} outer (entry side)"
        if dy - 40 <= x < dy + 40:
            return f"{name} inner"
        if dy + 40 <= x < hi:
            return f"{name} outer (exit side)"
    if N1[1] <= x < N2[0]:
        return "linker"
    return "before N1" if x < N1[0] else "after N2"


def nearest_landmark(x: float) -> tuple[str, float]:
    names = list(LANDMARKS)
    d = np.array([x - LANDMARKS[k] for k in names])
    j = int(np.argmin(np.abs(d)))
    return names[j], float(d[j])


# ----------------------------------------------------------------- data
def excellent_set() -> set[tuple[int, str]]:
    if TETHER_SET:
        return _user_set()
    c = pd.read_csv(p.DEST / "excellent_tethers" / "excellent_tethers_criteria.csv")
    c = c[c.selected.astype(str).str.lower().isin(("true", "1", "y"))]
    return {(int(d), str(b)) for d, b in zip(c.data_id, c.bead)}


def spike_mask(y: np.ndarray) -> np.ndarray:
    med = median_filter(y, size=2 * SPIKE_HALF + 1, mode="nearest")
    r = y - med
    s = 1.4826 * np.median(np.abs(r - np.median(r))) or 1.0
    return np.abs(r) <= SPIKE_NSIG * s


def collect() -> dict:
    exc = excellent_set()
    titr, _ = p.titration_files()
    plateaus, beads = [], []
    for d in titr:
        run = p.load(d)
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{d:03d}.h5")
        t_abs = run["time_s"]
        t = t_abs - t_abs[0]
        for roi in run["rois"]:
            bead = f"R{roi}"
            q = trk.get(roi, {})
            rec = {"data_id": d, "bead": bead, "tracking": q.get("status", "unknown"),
                   "level_jump_um": q.get("level_jump_um", np.nan), "excellent": (d, bead) in exc}
            beads.append(rec)
            if rec["tracking"] not in ("tracked", "lost") or not steps:
                rec["reason"] = f"tracking {rec['tracking']}"
                continue
            t_loss = (q["t_loss_s"] - t_abs[0]) if rec["tracking"] == "lost" else np.inf
            z_raw = run["z_um"][roi] * 1000.0
            z_s = F.smoothed_nm(run["z_um"][roi])
            run["z_um_raw"] = run["z_um"][roi]
            ext_s, note = F.rest_drift(run, steps, z_s)
            drift = z_s - ext_s
            ext_raw = z_raw - drift
            pl1 = steps[0]["plateau"]
            if abs(steps[0]["force_pn"] - ZERO_F) > 0.5:
                rec["reason"] = "first step not 10 pN"
                continue
            zero = F.max_extension(t[pl1], ext_s[pl1])
            rec.update(zero_nm=zero, drift=note)
            if not np.isfinite(zero):
                rec["reason"] = "no 10 pN zero"
                continue
            rec["reason"] = ""
            for st in steps:
                pl = st["plateau"]
                if t[pl[-1]] >= t_loss:
                    continue
                y = ext_raw[pl] - zero
                tt = t[pl]
                ok = np.isfinite(y)
                if ok.sum() < 50:
                    continue
                y, tt = y[ok], tt[ok]
                keep = spike_mask(y)
                up = st["up"]
                plateaus.append({
                    "data_id": d, "bead": bead, "step": st["step"], "force_pn": round(st["force_pn"], 1),
                    "pos_mm": st["pos_mm"], "excellent": (d, bead) in exc, "tracking": rec["tracking"],
                    "t": (tt[keep] - tt[0]).astype(float), "y": y[keep].astype(float),
                    "t_file_s": float(tt[0]),
                    "n_spikes": int((~keep).sum()), "hold_s": float(tt[-1] - tt[0]),
                    "fps": float(1.0 / np.median(np.diff(tt))),
                    "t_up_s": float(t[pl[0]] - t[up[0]]) if up.size else np.nan,
                    "zero_nm": zero, "drift": note})
        print(f"data_{d:03d}: {len(run['rois'])} ROIs, "
              f"{sum(1 for q in plateaus if q['data_id'] == d)} plateaus", flush=True)
    return {"plateaus": plateaus, "beads": pd.DataFrame(beads)}


def load_cache(recollect: bool = False) -> dict:
    if CACHE.exists() and not recollect:
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)
    OUT.mkdir(parents=True, exist_ok=True)
    data = collect()
    with open(CACHE, "wb") as fh:
        pickle.dump(data, fh)
    return data


# ----------------------------------------------------------------- detectors
def sigma_diff(y: np.ndarray) -> float:
    d = np.diff(y)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)) or 1.0


def op_segment(y: np.ndarray, sigma: float, lam: float, mmin: int = MMIN) -> list[tuple[int, int]]:
    """Exact optimal partitioning (PELT pruned): min SSE/sigma^2 + lam ln(n) per
    change point, every level >= mmin samples. Returns [(a, b)) in order."""
    n = y.size
    if n < 2 * mmin:
        return [(0, n)]
    beta = lam * np.log(n)
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y * y)])
    s2 = sigma * sigma
    Fv = np.full(n + 1, np.inf)
    Fv[0] = -beta
    last = np.zeros(n + 1, int)
    cand = np.array([0], int)
    for t in range(mmin, n + 1):
        s_new = t - mmin
        if s_new > 0 and np.isfinite(Fv[s_new]):
            cand = np.append(cand, s_new)
        if cand.size == 0:
            continue
        m = t - cand
        sse = (cs2[t] - cs2[cand]) - (cs[t] - cs[cand]) ** 2 / m
        tot = Fv[cand] + sse / s2 + beta
        j = int(np.argmin(tot))
        Fv[t], last[t] = tot[j], cand[j]
        cand = cand[(Fv[cand] + sse / s2) <= Fv[t]]
    out, t = [], n
    while t > 0:
        out.append((int(last[t]), int(t)))
        t = last[t]
    return out[::-1]


def bs_segment(y_nm: np.ndarray, f: float) -> list[tuple[int, int]]:
    """dwell_steps_20261009 exactly as written: BIC binary segmentation on Luz (model A)."""
    luz = (y_nm - handle_nm(f)[0]) / g_a(f)[0]
    sig = float(np.median(np.abs(np.diff(luz))) / 0.6745 / np.sqrt(2)) or 1.0
    return [(int(a), int(b)) for a, b in DS.merge(DS.segment(luz, sig), luz, DS.MERGE_BP)]


def segment(y, sigma, det, f):
    if det == "BS":
        return bs_segment(y, f)
    return op_segment(y, sigma, float(det))


def seg_table(t: np.ndarray, y: np.ndarray, segs, sigma: float, long_s: float = LONG_S) -> pd.DataFrame:
    n = y.size
    rows = []
    for k, (a, b) in enumerate(segs):
        t0 = t[a] if a == 0 else 0.5 * (t[a - 1] + t[a])
        t1 = t[b - 1] if b == n else 0.5 * (t[b - 1] + t[b])
        if a == 0 or b == n:                      # censored: add half a frame at the open end
            fr = float(np.median(np.diff(t))) if n > 1 else 0.0
            t0 = t0 - (0.5 * fr if a == 0 else 0.0)
            t1 = t1 + (0.5 * fr if b == n else 0.0)
        rows.append({"k": k, "a": a, "b": b, "t0": t0, "t1": t1, "duration_s": t1 - t0, "n_frames": b - a,
                     "level_nm": float(np.mean(y[a:b])), "sem_nm": sigma / np.sqrt(b - a),
                     "cens_start": a == 0, "cens_end": b == n})
    st = pd.DataFrame(rows)
    st["long"] = st.duration_s >= long_s
    return st


def classify(st: pd.DataFrame) -> pd.DataFrame:
    st = st.copy()
    st["category"] = np.where(st.long, "long", "outside")
    st["prev_long"] = -1
    st["next_long"] = -1
    st["returns"] = False
    longs = st.index[st.long].to_numpy()
    if longs.size < 2:
        return st
    cat, pl, nl, ret = st.category.to_numpy().copy(), st.prev_long.to_numpy().copy(), \
        st.next_long.to_numpy().copy(), st.returns.to_numpy().copy()
    lev, sem = st.level_nm.to_numpy(), st.sem_nm.to_numpy()
    for j in np.flatnonzero(~st.long.to_numpy()):
        pp, qq = longs[longs < j], longs[longs > j]
        if not pp.size or not qq.size:
            continue
        i1, i2 = int(pp[-1]), int(qq[0])
        l1, l2 = lev[i1], lev[i2]
        tol = max(MOVE_MIN_NM, 3.0 * np.hypot(sem[i1], sem[i2]))
        x = lev[j]
        same = abs(l2 - l1) <= tol
        if same or x > max(l1, l2) or x < min(l1, l2):
            c = "excursion_fwd" if x > (l1 if same else max(l1, l2)) else "excursion_back"
        else:
            c = "intermediate_fwd" if l2 > l1 else "intermediate_back"
        cat[j], pl[j], nl[j], ret[j] = c, i1, i2, same
    st["category"], st["prev_long"], st["next_long"], st["returns"] = cat, pl, nl, ret
    return st


def detect(t, y, sigma, det, f, long_s=LONG_S) -> pd.DataFrame:
    return classify(seg_table(t, y, segment(y, sigma, det, f), sigma, long_s))


def counts(st: pd.DataFrame) -> dict:
    c = st.category.value_counts()
    out = {k: int(c.get(k, 0)) for k in CATS}
    out["between"] = sum(out[k] for k in CATS)
    out["brief_all"] = int((~st.long).sum())
    out["n_long"] = int(st.long.sum())
    lg = st[st.long].level_nm.to_numpy()
    out["n_fwd_moves"] = int(np.sum(np.diff(lg) > MOVE_MIN_NM)) if lg.size > 1 else 0
    out["excursion_returns"] = int((st.category.str.startswith("excursion") & st.returns).sum())
    return out


# ----------------------------------------------------------------- nulls
def skeleton(st: pd.DataFrame, n: int) -> np.ndarray:
    """Long-stall skeleton: every brief level replaced by the nearer long stall in time."""
    longs = st[st.long]
    if not len(longs):
        return np.full(n, float(np.average(st.level_nm, weights=st.n_frames)))
    mids = 0.5 * (longs.a + longs.b).to_numpy()
    lv = longs.level_nm.to_numpy()
    sk = np.empty(n)
    for r in st.itertuples():
        if r.long:
            sk[r.a:r.b] = r.level_nm
        else:
            for i in range(r.a, r.b):                       # frame-wise nearest long stall
                sk[i] = lv[int(np.argmin(np.abs(mids - i)))]
    return sk


def creep_base(st: pd.DataFrame, t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """No discrete brief level anywhere: each long stall its own line (creep),
    consecutive long stalls joined by a straight ramp over the interval the
    brief levels occupied, flanks before / after as lines."""
    n = y.size
    base = np.full(n, np.nan)
    longs = list(st[st.long].itertuples())

    def line(a_, b_):
        if b_ - a_ >= 3:
            return np.polyval(np.polyfit(t[a_:b_], y[a_:b_], 1), t[a_:b_])
        return np.full(b_ - a_, float(np.mean(y[a_:b_])))
    if not longs:
        return line(0, n)
    for r in longs:
        base[r.a:r.b] = line(r.a, r.b)
    for r1, r2 in zip(longs[:-1], longs[1:]):
        if r2.a > r1.b:
            base[r1.b:r2.a] = np.interp(t[r1.b:r2.a], [t[r1.b - 1], t[r2.a]], [base[r1.b - 1], base[r2.a]])
    if longs[0].a > 0:
        base[:longs[0].a] = line(0, longs[0].a)
    if longs[-1].b < n:
        base[longs[-1].b:] = line(longs[-1].b, n)
    return base


def ar1_noise(n: int, sd: float, rho: float, rng) -> np.ndarray:
    rho = float(np.clip(rho, -0.5, 0.8))
    e = rng.normal(0.0, sd * np.sqrt(1 - rho * rho), n)
    x = np.empty(n)
    x[0] = rng.normal(0.0, sd)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + e[i]
    return x


def block_boot(res: np.ndarray, n: int, rng, block: int = BLOCK) -> np.ndarray:
    m = res.size
    k = int(np.ceil(n / block))
    starts = rng.integers(0, m, k)
    idx = (starts[:, None] + np.arange(block)[None, :]) % m
    return res[idx.ravel()[:n]]


def residuals(st: pd.DataFrame, y: np.ndarray) -> np.ndarray:
    fit = np.concatenate([np.full(r.b - r.a, r.level_nm) for r in st.itertuples()])
    return y - fit


# ----------------------------------------------------------------- QC (stage prep)
_P: list = []


def _init(plats):
    global _P
    _P = plats


def qc_one(i: int) -> dict:
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    sig = sigma_diff(y)
    r9 = y - median_filter(y, 9, mode="nearest")
    _, cnt = np.unique(np.round(y, 3), return_counts=True)
    st = detect(t, y, sig, LAM_REF, f)
    res = residuals(st, y)
    rho = float(np.corrcoef(res[:-1], res[1:])[0, 1]) if res.std() > 0 else 0.0
    lv = abs_bp(st.level_nm.to_numpy(), f, "A")
    rec = {"i": i, "data_id": q["data_id"], "bead": q["bead"], "step": q["step"], "force_pn": f,
           "excellent": q["excellent"], "tracking": q["tracking"], "fps": round(q["fps"], 1),
           "n_frames": y.size, "hold_s": round(q["hold_s"], 2), "n_spikes": q["n_spikes"],
           "sigma_nm": round(sig, 2), "sigma_res_nm": round(float(res.std()), 2), "rho1": round(rho, 3),
           "acf": [round(float(np.corrcoef(res[:-k], res[k:])[0, 1]), 3) if res.std() > 0 else 0.0
                   for k in range(1, 11)],
           "flicker": round(float(np.mean(np.abs(r9) > 4 * sig)), 4),
           "repeat": round(float(cnt.max() / y.size), 4),
           "lo_bp_A": round(float(lv.min()), 1), "hi_bp_A": round(float(lv.max()), 1)}
    reason = []
    if not SIG_MIN_NM <= sig <= SIG_MAX_NM:
        reason.append("sigma")
    if rec["flicker"] > FLICKER_MAX:
        reason.append("flicker")
    if rec["repeat"] > REPEAT_MAX:
        reason.append("clipped")
    if lv.min() < GUARD_BP[0] or lv.max() > GUARD_BP[1]:
        reason.append("LUT level")
    rec["gate_reason"] = "+".join(reason)
    rec["passes"] = not reason
    return rec


# ----------------------------------------------------------------- scan (stage scan)
DETECTORS = tuple(LAMBDAS) + ("BS",)


def det_name(det) -> str:
    return "BS (dwell_steps)" if det == "BS" else f"OP lambda={float(det):g}"


def scan_one(i: int) -> tuple[list[dict], list[dict]]:
    """Real trace and N1 / N2 / N3 nulls, every detector. Returns count rows and
    event rows (brief between-long events with amplitudes) for the nulls."""
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    rng = np.random.default_rng(1000 + i)
    sig = sigma_diff(y)
    st = detect(t, y, sig, LAM_REF, f)
    res = residuals(st, y)
    sres = float(res.std())
    rho = float(np.corrcoef(res[:-1], res[1:])[0, 1]) if sres > 0 else 0.0
    sk = skeleton(st, y.size)
    cr = creep_base(st, t, y)
    traces = [("real", 0, y)]
    for rep in range(N_NULL):
        traces.append(("N1 white", rep, sk + ar1_noise(y.size, sres, rho, rng)))
        traces.append(("N2 boot", rep, sk + block_boot(res, y.size, rng)))
        traces.append(("N3 creep", rep, cr + block_boot(res, y.size, rng)))
    crow, erow = [], []
    for name, rep, ys in traces:
        s_ = sigma_diff(ys)
        for det in DETECTORS:
            sst = detect(t, ys, s_, det, f)
            crow.append({"i": i, "force_pn": f, "excellent": q["excellent"], "null": name, "rep": rep,
                         "detector": det_name(det), "plateau_s": q["hold_s"], **counts(sst)})
            if name == "real":
                continue
            lev = sst.level_nm.to_numpy()
            for r in sst[sst.category.isin(CATS)].itertuples():
                erow.append({"i": i, "force_pn": f, "null": name, "rep": rep, "detector": det_name(det),
                             "category": r.category, "returns": r.returns, "duration_s": r.duration_s,
                             "n_frames": r.n_frames, "amp_prev_nm": r.level_nm - lev[r.prev_long],
                             "gap_nm": lev[r.next_long] - lev[r.prev_long]})
    return crow, erow


# ----------------------------------------------------------------- injection (stage inject)
# The camera clock is IRREGULAR: inside one plateau the frame spacing runs from
# 0.01 to ~0.19 s (median 0.01 / 0.02 / 0.03-0.04 s, mean rate well below the
# median), so a brief level can fall into a gap. Injections therefore use REAL
# plateau time grids, grouped by their median frame rate, and durations in s.
INJ_CLASSES = {"100 fps": (90.0, 1e9), "50 fps": (40.0, 90.0), "<=33 fps": (0.0, 40.0)}
INJ_SIG = (25.0, 36.0, 50.0)
INJ_DURS = (0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
INJ_AMPS = (20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0)
INJ_REPS = 40
N_GRIDS = 150


def inject_one(args) -> list[dict]:
    cls, grids, sig, det, seed = args
    rng = np.random.default_rng(seed)
    rows = []
    for dur in INJ_DURS:
        for amp in INJ_AMPS:
            hit = {"fwd": 0, "exc": 0}
            err = {"fwd": [], "exc": []}
            nfr = []
            for rep in range(INJ_REPS):
                t = grids[int(rng.integers(len(grids)))]
                n = t.size
                t0 = rng.uniform(2.5, min(6.0, t[-1] - dur - 2.0))
                inj = np.flatnonzero((t >= t0) & (t < t0 + dur))
                nfr.append(inj.size)
                if inj.size == 0:
                    continue                         # the event fell into a frame gap: missed
                noise = rng.normal(0.0, sig, n)
                for kind in ("fwd", "exc"):
                    base = np.zeros(n)
                    base[inj] = amp
                    if kind == "fwd":
                        base[inj[-1] + 1:] = 2 * amp
                    ys = base + noise
                    st = detect(t, ys, sigma_diff(ys), det, 30.0)
                    want = "intermediate_fwd" if kind == "fwd" else "excursion_fwd"
                    for r in st[st.category == want].itertuples():
                        ov = np.sum((inj >= r.a) & (inj < r.b))
                        if ov >= 0.5 * inj.size and abs(r.level_nm - amp) <= amp / 2:
                            hit[kind] += 1
                            err[kind].append((r.level_nm - amp, r.duration_s - dur))
                            break
            for kind in ("fwd", "exc"):
                e = np.array(err[kind]) if err[kind] else np.full((1, 2), np.nan)
                rows.append({"fps_class": cls, "sigma_nm": sig, "detector": det_name(det), "kind": kind,
                             "duration_s": dur, "amp_nm": amp, "mean_frames_in_event": round(float(np.mean(nfr)), 2),
                             "frac_event_in_gap": round(float(np.mean(np.array(nfr) == 0)), 3),
                             "recall": hit[kind] / INJ_REPS, "amp_bias_nm": float(np.nanmedian(e[:, 0])),
                             "dur_bias_s": float(np.nanmedian(e[:, 1]))})
    return rows


def ramp_vs_steps(st, i1, i2, y, sigma, t, pad=25, return_fit=False):
    """BIC(ramp) - BIC(steps) over L1 tail .. L2 head (units sigma^2); > 0 favours
    the discrete intermediate level(s) over a continuous ramp, LINEAR IN TIME
    (irregular frame clock), whose two breakpoints are free (a single step is the
    limit of a ramp, so dBIC > 0 also beats a direct jump)."""
    a = max(int(st.a[i1]), int(st.b[i1]) - pad)
    b = min(int(st.b[i2]), int(st.a[i2]) + pad)
    w = y[a:b]
    n = w.size
    fit = np.concatenate([np.full(int(min(r.b, b) - max(r.a, a)), r.level_nm)
                          for r in st.loc[i1:i2].itertuples() if min(r.b, b) > max(r.a, a)])
    k_lv = i2 - i1 + 1
    p_steps = k_lv + (k_lv - 1)
    sse_s = float(np.sum((w - fit) ** 2))
    x = np.asarray(t[a:b], float)
    best, best_fit = np.inf, None
    for c1 in range(0, n - 1):
        c2 = np.arange(c1 + 1, n)[:, None]
        h = np.clip((x[None, :] - x[c1]) / (x[c2] - x[c1]), 0.0, 1.0)
        u = 1.0 - h
        s11, s12, s22 = (u * u).sum(1), (u * h).sum(1), (h * h).sum(1)
        b1, b2 = (u * w).sum(1), (h * w).sum(1)
        det = s11 * s22 - s12 * s12
        ok = det > 1e-9
        l1 = np.where(ok, (s22 * b1 - s12 * b2) / np.where(ok, det, 1), np.nan)
        l2 = np.where(ok, (s11 * b2 - s12 * b1) / np.where(ok, det, 1), np.nan)
        sse = ((w[None, :] - (l1[:, None] * u + l2[:, None] * h)) ** 2).sum(1)
        if np.isfinite(sse).any() and float(np.nanmin(sse)) < best:
            j = int(np.nanargmin(sse))
            best = float(sse[j])
            if return_fit:
                best_fit = l1[j] * u[j] + l2[j] * h[j]
    dbic = (best / sigma ** 2 + 4 * np.log(n)) - (sse_s / sigma ** 2 + p_steps * np.log(n))
    return (dbic, (a, b, best_fit)) if return_fit else dbic


def events_one(args) -> tuple[list[dict], list[dict], list[dict], dict]:
    i, lam, long_s = args
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    ga, gb = gains(f)
    sig = sigma_diff(y)
    st = detect(t, y, sig, lam, f, long_s)
    lev = st.level_nm.to_numpy()
    base = {"i": i, "data_id": q["data_id"], "bead": q["bead"], "step": q["step"], "force_pn": f,
            "excellent": q["excellent"], "fps": round(q["fps"], 1), "sigma_nm": round(sig, 1),
            "t_file_s": round(q["t_file_s"], 2)}
    ev, lg, tr = [], [], []
    for r in st[st.long].itertuples():
        lg.append({**base, "k": r.k, "t_start_s": round(r.t0, 3), "duration_s": round(r.duration_s, 3),
                   "censored": bool(r.cens_start or r.cens_end), "level_nm": round(r.level_nm, 1),
                   "abs_bp_A": round(float(abs_bp(r.level_nm, f, "A")[0]), 1),
                   "abs_bp_B": round(float(abs_bp(r.level_nm, f, "B")[0]), 1)})
    longs = st.index[st.long].to_numpy()
    dbic_cache = {}
    for i1, i2 in zip(longs[:-1], longs[1:]):
        l1, l2 = lev[i1], lev[i2]
        tol = max(MOVE_MIN_NM, 3.0 * np.hypot(st.sem_nm[i1], st.sem_nm[i2]))
        mid = st.loc[i1 + 1:i2 - 1]
        direction = "same" if abs(l2 - l1) <= tol else ("fwd" if l2 > l1 else "back")
        n_int = int(mid.category.str.startswith("intermediate").sum())
        db = np.nan
        if direction != "same" and n_int:
            db = float(ramp_vs_steps(st, int(i1), int(i2), y, sig, t))
        dbic_cache[(int(i1), int(i2))] = db
        tr.append({**base, "k1": int(i1), "k2": int(i2), "l1_nm": round(l1, 1), "l2_nm": round(l2, 1),
                   "gap_nm": round(l2 - l1, 1), "gap_bp_A": round((l2 - l1) / ga, 1),
                   "gap_bp_B": round((l2 - l1) / gb, 1), "direction": direction,
                   "l1_s": round(st.duration_s[i1], 3), "l2_s": round(st.duration_s[i2], 3),
                   "t_l1_end_s": round(st.t1[i1], 3), "t_l2_start_s": round(st.t0[i2], 3),
                   "n_brief": len(mid), "n_intermediate": n_int,
                   "n_excursion": int(mid.category.str.startswith("excursion").sum()),
                   "l1_abs_bp_A": round(float(abs_bp(l1, f, "A")[0]), 1),
                   "l2_abs_bp_A": round(float(abs_bp(l2, f, "A")[0]), 1),
                   "l1_abs_bp_B": round(float(abs_bp(l1, f, "B")[0]), 1),
                   "l2_abs_bp_B": round(float(abs_bp(l2, f, "B")[0]), 1),
                   "dbic_ramp_minus_steps": round(db, 1) if np.isfinite(db) else np.nan})
    for r in st[~st.long].itertuples():
        row = {**base, "k": r.k, "category": r.category, "returns": bool(r.returns),
               "t_start_s": round(r.t0, 3), "duration_s": round(r.duration_s, 3), "n_frames": r.n_frames,
               "censored": bool(r.cens_start or r.cens_end), "level_nm": round(r.level_nm, 1),
               "sem_nm": round(r.sem_nm, 1),
               "jump_in_nm": round(r.level_nm - lev[r.k - 1], 1) if r.k > 0 else np.nan,
               "jump_out_nm": round(lev[r.k + 1] - r.level_nm, 1) if r.k < len(st) - 1 else np.nan}
        for mdl, g in (("A", ga), ("B", gb)):
            x = float(abs_bp(r.level_nm, f, mdl)[0])
            nm_, dd = nearest_landmark(x)
            row.update({f"abs_bp_{mdl}": round(x, 1), f"region_{mdl}": turn_region(x),
                        f"nearest_landmark_{mdl}": nm_, f"to_landmark_bp_{mdl}": round(dd, 1),
                        f"jump_in_bp_{mdl}": round(row["jump_in_nm"] / g, 1) if np.isfinite(row["jump_in_nm"]) else np.nan,
                        f"jump_out_bp_{mdl}": round(row["jump_out_nm"] / g, 1) if np.isfinite(row["jump_out_nm"]) else np.nan})
        if r.prev_long >= 0:
            l1, l2 = lev[r.prev_long], lev[r.next_long]
            d1, d2, gap = r.level_nm - l1, l2 - r.level_nm, l2 - l1
            n_in_gap = int(((st.k > r.prev_long) & (st.k < r.next_long)).sum())
            row.update({"prev_long_nm": round(l1, 1), "next_long_nm": round(l2, 1),
                        "prev_long_s": round(st.duration_s[r.prev_long], 3),
                        "next_long_s": round(st.duration_s[r.next_long], 3),
                        "prev_long_end_s": round(st.t1[r.prev_long], 3),
                        "next_long_start_s": round(st.t0[r.next_long], 3),
                        "n_brief_in_gap": n_in_gap,
                        "d1_nm": round(d1, 1), "d2_nm": round(d2, 1), "gap_nm": round(gap, 1),
                        "u": round(d1 / gap, 3) if abs(gap) > 1e-9 and not r.returns else np.nan,
                        "amp_nm": round(r.level_nm - (0.5 * (l1 + l2) if r.returns else l1), 1),
                        "dbic_ramp_minus_steps": dbic_cache.get((int(r.prev_long), int(r.next_long)), np.nan)})
            for mdl, g in (("A", ga), ("B", gb)):
                row.update({f"d1_bp_{mdl}": round(d1 / g, 1), f"d2_bp_{mdl}": round(d2 / g, 1),
                            f"gap_bp_{mdl}": round(gap / g, 1), f"amp_bp_{mdl}": round(row["amp_nm"] / g, 1),
                            f"prev_long_abs_bp_{mdl}": round(float(abs_bp(l1, f, mdl)[0]), 1),
                            f"next_long_abs_bp_{mdl}": round(float(abs_bp(l2, f, mdl)[0]), 1)})
            # both flanking main stalls at the start (< 100 bp, model A): the fork
            # has not reached N1 (entry 314 bp), so the event cannot be a turn transition
            row["at_closed_fork"] = bool(max(row["prev_long_abs_bp_A"], row["next_long_abs_bp_A"]) < 100)
        ev.append(row)
    cnt = {**base, "plateau_s": q["hold_s"], **counts(st),
           "max_long_bp_A": round(float(abs_bp(st[st.long].level_nm.max(), f, "A")[0]), 1)
           if st.long.any() else np.nan}
    return ev, lg, tr, cnt


# ----------------------------------------------------------------- geometry nulls (stage geom)
def sim_window(t_grid, sigma, levels, durs, rng, lam, res=None):
    """Synthetic window on a REAL time grid (a random stretch of the plateau's
    own clock): consecutive levels (nm) held for `durs` (s); noise = block
    bootstrap of the plateau's residuals (or Gaussian)."""
    total = float(np.sum(durs))
    t0 = rng.uniform(0.0, max(t_grid[-1] - total, 0.0))
    w = (t_grid >= t0) & (t_grid < t0 + total)
    t = t_grid[w] - t0
    n = t.size
    edges = np.cumsum(durs)[:-1]
    base = np.asarray(levels, float)[np.searchsorted(edges, t, side="right")]
    noise = block_boot(res, n, rng) if res is not None else rng.normal(0.0, sigma, n)
    ys = base + noise
    s_ = sigma_diff(ys)
    st = detect(t, ys, s_, lam, 30.0)
    st.attrs["y"], st.attrs["sigma"] = ys, s_
    return st


def geom_one(args) -> list[dict]:
    """Detector-aware uniform null for one real forward transition with intermediates:
    one intermediate at a uniform fraction u of the gap, duration drawn from the
    observed brief durations, own grid / residual noise / long-stall durations."""
    (i, gap_nm, l1_s, l2_s, durs, lam, l1_abs_A, l1_abs_B, f, seed, mode) = args
    q = _P[i]
    rng = np.random.default_rng(seed)
    y, t = q["y"], q["t"]
    sig = sigma_diff(y)
    st0 = detect(t, y, sig, LAM_REF, f)
    res = residuals(st0, y)
    ga, gb = gains(f)
    rows = []
    for rep in range(N_GEOM):
        tau = float(rng.choice(durs))
        if mode is None:                    # uniform position inside the gap
            u = rng.uniform(0.0, 1.0)
        else:                               # H1 'ladder': intermediate k x period bp after L1
            _, per, mdl = mode
            g = ga if mdl == "A" else gb
            kmax = int(np.floor((gap_nm / g - 1e-6) / per))
            if kmax < 1:
                continue
            u = rng.integers(1, kmax + 1) * per * g / gap_nm
            if u >= 1.0:
                continue
        st = sim_window(t, sig, [0.0, u * gap_nm, gap_nm],
                        [float(np.clip(l1_s, 1.05, 3.0)), tau, float(np.clip(l2_s, 1.05, 3.0))], rng, lam, res)
        for r in st[st.category == "intermediate_fwd"].itertuples():
            l1, l2 = st.level_nm[int(r.prev_long)], st.level_nm[int(r.next_long)]
            d1 = r.level_nm - l1
            rows.append({"i": i, "rep": rep, "mode": "uniform" if mode is None else f"ladder {mode[1]:g} bp {mode[2]}",
                         "u_true": u, "tau_s": tau, "force_pn": f,
                         "u": d1 / (l2 - l1) if l2 != l1 else np.nan,
                         "d1_bp_A": d1 / ga, "d1_bp_B": d1 / gb,
                         "d2_bp_A": (l2 - r.level_nm) / ga, "d2_bp_B": (l2 - r.level_nm) / gb,
                         "abs_bp_A": l1_abs_A + d1 / ga, "abs_bp_B": l1_abs_B + d1 / gb,
                         "duration_s": r.duration_s})
    return rows


# ----------------------------------------------------------------- driver helpers
def pool_map(fn, items, plats, procs=9, chunksize=4):
    ctx = mp.get_context("fork")
    with ctx.Pool(procs, initializer=_init, initargs=(plats,)) as pool:
        return pool.map(fn, items, chunksize=chunksize)


# ================================================================= stages
SETS = {"excellent": lambda d: d.excellent.astype(bool), "all_clean": lambda d: np.ones(len(d), bool)}
SET_COL = {"excellent": "#E64B35", "all_clean": "#3C5488"}


def stage_prep(plats) -> pd.DataFrame:
    rows = pool_map(qc_one, list(range(len(plats))), plats)
    qc = pd.DataFrame(rows)
    acf = np.vstack(qc.pop("acf").to_numpy())
    for k in range(acf.shape[1]):
        qc[f"acf{k + 1}"] = acf[:, k]
    qc.to_csv(OUT / "plateaus_qc.csv", index=False)
    print(f"prep: {len(qc)} plateaus, {int(qc.passes.sum())} pass "
          f"({int(qc[qc.passes].groupby(['data_id', 'bead']).ngroups)} tethers; excellent "
          f"{int(qc[qc.passes & qc.excellent].groupby(['data_id', 'bead']).ngroups)} tethers, "
          f"{int((qc.passes & qc.excellent).sum())} plateaus)", flush=True)
    return qc


def stage_scan(plats, qc) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    idx = qc[qc.passes].i.tolist()
    res = pool_map(scan_one, idx, plats, chunksize=2)
    cnt = pd.DataFrame([r for c, _ in res for r in c])
    evn = pd.DataFrame([r for _, e in res for r in e])
    cnt.to_csv(OUT / "detector_null_scan_counts.csv.gz", index=False)
    evn.to_csv(OUT / "detector_null_scan_events.csv.gz", index=False)
    lam = choose_lambda(cnt)
    return cnt, evn, lam


def scan_summary(cnt: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (det, nul, hi), g in cnt.assign(hi=np.where(cnt.force_pn >= 18, ">=18 pN", "10 pN")) \
            .groupby(["detector", "null", "hi"]):
        reps = g.rep.nunique()
        nplat = g.i.nunique()
        rows.append({"detector": det, "null": nul, "forces": hi, "n_plateaus": nplat,
                     **{f"{c}_per_plateau": round(g[c].sum() / reps / nplat, 4)
                        for c in ("between",) + CATS + ("brief_all",)},
                     "n_long_per_plateau": round(g.n_long.sum() / reps / nplat, 3)})
    return pd.DataFrame(rows)


def choose_lambda(cnt: pd.DataFrame) -> float:
    s = scan_summary(cnt)
    s.to_csv(OUT / "detector_null_scan_summary.csv", index=False)
    hi = s[s.forces == ">=18 pN"].set_index(["detector", "null"]).between_per_plateau
    # lambda is set against the NOISE nulls (N1, N2): is a brief level a noise
    # artefact? N3 (continuous ramp over the same interval) asks a different
    # question - discrete level or finite-speed move - and is answered per
    # transition by BIC(ramp) - BIC(steps) in stage 'disc', not by counts.
    pick = None
    for lam in LAMBDAS:
        d = det_name(lam)
        real = hi[(d, "real")]
        worst = max(hi[(d, n)] for n in ("N1 white", "N2 boot"))
        print(f"  {d}: real between-long {real:.4f}/plateau, worst noise null {worst:.4f} "
              f"({worst / max(real, 1e-9):.1%}); N3 creep {hi[(d, 'N3 creep')]:.4f}", flush=True)
        if pick is None and worst <= FP_TARGET_FRAC * real:
            pick = lam
    d = det_name("BS")
    print(f"  {d}: real {hi[(d, 'real')]:.3f}, nulls "
          + ", ".join(f"{n} {hi[(d, n)]:.4f}" for n in ("N1 white", "N2 boot", "N3 creep")), flush=True)
    return float(pick if pick is not None else LAMBDAS[-1])


def stage_inject(lam: float, plats, qc) -> pd.DataFrame:
    g = qc[qc.passes & (qc.force_pn >= 18) & (qc.hold_s >= 9.5)]
    rng = np.random.default_rng(3)
    grids = {}
    for cls, (lo, hi) in INJ_CLASSES.items():
        idx = g[(g.fps >= lo) & (g.fps < hi)].i.to_numpy()
        idx = rng.choice(idx, min(N_GRIDS, idx.size), replace=False)
        grids[cls] = [plats[i]["t"][plats[i]["t"] <= 10.0] for i in idx]
    jobs = [(cls, grids[cls], sig, det, 7 + 31 * k) for k, (cls, sig, det) in enumerate(
        (c_, s_, d_) for c_ in INJ_CLASSES for s_ in INJ_SIG for d_ in (lam, "BS"))]
    res = pool_map(inject_one, jobs, [], chunksize=1)
    inj = pd.DataFrame([r for rr in res for r in rr])
    inj.to_csv(OUT / "detector_injection_recall.csv", index=False)
    return inj


def detection_limits(inj: pd.DataFrame) -> pd.DataFrame:
    """Shortest brief dwell (< LONG_S; longer ones are long stalls by definition)
    from which recall stays >= 0.8, per amplitude, frame rate, noise, detector."""
    rows = []
    inj = inj[inj.duration_s < 0.95 * LONG_S]
    for (fps, sig, det, kind), g in inj.groupby(["fps_class", "sigma_nm", "detector", "kind"]):
        for amp, ga_ in g.groupby("amp_nm"):
            mins = np.nan
            rr = ga_.sort_values("duration_s")
            for j in range(len(rr)):
                if (rr.recall.iloc[j:] >= 0.8).all():
                    mins = rr.duration_s.iloc[j]
                    break
            rows.append({"fps_class": fps, "sigma_nm": sig, "detector": det, "kind": kind, "amp_nm": amp,
                         "amp_bp_A_30pN": round(amp / gains(30.0)[0], 1),
                         "amp_bp_B_30pN": round(amp / gains(30.0)[1], 1),
                         "shortest_dwell_recall80_s": mins,
                         "recall_at_0.1s": float(np.interp(0.1, rr.duration_s, rr.recall)),
                         "recall_at_0.3s": float(np.interp(0.3, rr.duration_s, rr.recall))})
    lim = pd.DataFrame(rows)
    lim.to_csv(OUT / "detector_detection_limits.csv", index=False)
    return lim


def channel_flags(ev: pd.DataFrame, plats) -> pd.DataFrame:
    """Tracking check for every brief level between long stalls: robust z of the
    in-level median of X, Y and the LUT amplitude A against +-1 s around it.
    A Z event that moves X / Y / A by > 6 robust sigma is flagged as a tracking
    suspect. Large fast excursions (|amp| >= 150 nm in <= 5 frames, i.e. a
    200-550 bp opening and re-zipping within 0.1 s at constant force) are
    flagged too."""
    import h5py
    ev = ev.copy()
    for c in ("zX", "zY", "zA"):
        ev[c] = np.nan
    sel = ev.category.isin(CATS)
    for (d, bead, step), g in ev[sel].groupby(["data_id", "bead", "step"]):
        run = p.load(d)
        steps = p.segment_steps(run)
        pl = steps[step - 1]["plateau"]
        t = run["time_s"]
        with h5py.File(p.RAW / f"data_{d:03d}.h5") as f:
            ch = {c: np.asarray(f[f"Traces/{nm_}"][bead[1:]], float)
                  for c, nm_ in (("X", "X (um)"), ("Y", "Y (um)"), ("A", "A (a.u.)"))}
        for j, r in g.iterrows():
            ts = t[pl[0]] + r.t_start_s
            te = ts + r.duration_s
            w_in = (t >= ts) & (t < te)
            w_out = ((t >= ts - 1) & (t < ts)) | ((t >= te) & (t < te + 1))
            for c, v in ch.items():
                o = v[w_out]
                o = o[np.isfinite(o)]
                vi = v[w_in]
                vi = vi[np.isfinite(vi)]
                if o.size < 5 or vi.size < 1:
                    continue
                mad = 1.4826 * np.median(np.abs(o - np.median(o))) or 1e-12
                ev.loc[j, f"z{c}"] = round(float((np.median(vi) - np.median(o)) / mad), 1)
    zz = ev[["zX", "zY", "zA"]].abs().max(axis=1)
    ev["flag_tracking"] = sel & (zz > 6)
    ev["flag_large_fast"] = sel & ev.category.str.startswith("excursion") & \
        (ev.amp_nm.abs() >= 150) & (ev.n_frames <= 5)
    return ev


def stage_events(plats, qc, lam) -> dict:
    idx = qc[qc.passes].i.tolist()
    out = {}
    for long_s in (LONG_S,) + LONG_SENS:
        res = pool_map(events_one, [(i, lam, long_s) for i in idx], plats, chunksize=4)
        ev = pd.DataFrame([r for e, _, _, _ in res for r in e])
        lg = pd.DataFrame([r for _, l_, _, _ in res for r in l_])
        tr = pd.DataFrame([r for _, _, t_, _ in res for r in t_])
        pc = pd.DataFrame([c for _, _, _, c in res])
        if long_s == LONG_S:
            ev = channel_flags(ev, plats)
            ev.to_csv(OUT / "events_brief_all.csv", index=False)
            ev[ev.category.isin(CATS)].to_csv(OUT / "events_between_long.csv", index=False)
            lg.to_csv(OUT / "long_stalls.csv", index=False)
            tr.to_csv(OUT / "transitions_long_to_long.csv", index=False)
            pc.to_csv(OUT / "plateau_counts.csv", index=False)
        out[long_s] = (ev, lg, tr, pc)
        print(f"events (long >= {long_s:g} s): {len(ev)} brief levels, "
              f"{int(ev.category.isin(CATS).sum())} between long stalls, {len(lg)} long stalls", flush=True)
    return out


def long_threshold_table(out: dict) -> pd.DataFrame:
    rows = []
    for long_s, (ev, lg, tr, pc) in out.items():
        for sname, sel in SETS.items():
            for hi, fsel in (("10 pN", lambda d: d.force_pn == 10.0), (">=18 pN", lambda d: d.force_pn >= 18)):
                e = ev[sel(ev) & fsel(ev)]
                c = pc[sel(pc) & fsel(pc)]
                b = e[e.category.isin(CATS)]
                rows.append({"long_s": long_s, "set": sname, "forces": hi, "n_plateaus": len(c),
                             "long_per_plateau": round(c.n_long.mean(), 3),
                             "between_per_plateau": round(len(b) / max(len(c), 1), 4),
                             **{f"n_{k}": int((b.category == k).sum()) for k in CATS},
                             "median_brief_s": round(float(b.duration_s.median()), 3) if len(b) else np.nan})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- statistics
def comb_c(x, period):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.mean(np.cos(2 * np.pi * x / period))) if x.size else np.nan


def landmark_dist(x):
    x = np.asarray(x, float)
    return np.min(np.abs(x[:, None] - LM[None, :]), axis=1)


def ladder_table(ev, tr) -> pd.DataFrame:
    """One row per forward main-stall transition with >= 1 brief intermediate:
    the level ladder L1 -> I1 -> ... -> L2 as consecutive steps (nm, bp A, bp B)."""
    t = tr[(tr.direction == "fwd") & (tr.n_intermediate > 0)]
    rows = []
    for r in t.itertuples():
        e = ev[(ev.i == r.i) & (ev.prev_long_nm == r.l1_nm) & (ev.next_long_nm == r.l2_nm)].sort_values("t_start_s")
        lv = np.array([r.l1_nm] + e.level_nm.tolist() + [r.l2_nm])
        ga, gb = gains(r.force_pn)
        st_ = np.diff(lv)
        rows.append({"data_id": r.data_id, "bead": r.bead, "step": r.step, "force_pn": r.force_pn,
                     "excellent": r.excellent, "sigma_nm": r.sigma_nm, "fps": r.fps,
                     "L1_abs_bp_A": r.l1_abs_bp_A, "L2_abs_bp_A": r.l2_abs_bp_A,
                     "L1_abs_bp_B": r.l1_abs_bp_B, "L2_abs_bp_B": r.l2_abs_bp_B,
                     "gap_bp_A": r.gap_bp_A, "gap_bp_B": r.gap_bp_B, "n_brief": len(e),
                     "t_L1_end_s": r.t_l1_end_s, "t_L2_start_s": r.t_l2_start_s,
                     "steps_nm": " ".join(f"{x:.0f}" for x in st_),
                     "steps_bp_A": " ".join(f"{x / ga:.0f}" for x in st_),
                     "steps_bp_B": " ".join(f"{x / gb:.0f}" for x in st_),
                     "brief_durations_s": " ".join(f"{x:.2f}" for x in e.duration_s),
                     "brief_abs_bp_A": " ".join(f"{x:.0f}" for x in e.abs_bp_A),
                     "brief_abs_bp_B": " ".join(f"{x:.0f}" for x in e.abs_bp_B),
                     "dbic_ramp_minus_steps": r.dbic_ramp_minus_steps,
                     "discrete_supported": bool(r.dbic_ramp_minus_steps > DBIC_STRONG)})
    return pd.DataFrame(rows).sort_values("dbic_ramp_minus_steps", ascending=False)


def stage_geom(plats, ev, tr, lam) -> pd.DataFrame:
    """Detector-aware uniform null for every forward transition with >= 1 forward intermediate."""
    fi = ev[ev.category == "intermediate_fwd"]
    durs = fi.duration_s.to_numpy()
    jobs = []
    modes = [None] + [("ladder", per, mdl) for per in (40.0, 80.0) for mdl in ("A", "B")]
    for k, r in enumerate(tr[(tr.direction == "fwd") & (tr.n_intermediate > 0)].itertuples()):
        for m_, mode in enumerate(modes):
            jobs.append((r.i, r.gap_nm, r.l1_s, r.l2_s, durs, lam, r.l1_abs_bp_A, r.l1_abs_bp_B,
                         r.force_pn, 50_000 + 1000 * m_ + k, mode))
    res = pool_map(geom_one, jobs, plats, chunksize=2)
    un = pd.DataFrame([r for rr in res for r in rr])
    un["excellent"] = [plats[i]["excellent"] for i in un.i]
    un.to_csv(OUT / "geometry_null_uniform_detector_aware.csv.gz", index=False)
    return un


def geometry_tests(ev, un, evn, lam, rng_seed=11) -> pd.DataFrame:
    from scipy.stats import ks_2samp, mannwhitneyu
    rng = np.random.default_rng(rng_seed)
    rows = []
    det = det_name(lam)
    fp = evn[(evn.detector == det) & (evn.force_pn >= 18)]
    lad = un[un["mode"] != "uniform"]
    un = un[un["mode"] == "uniform"]
    for sname, sel in SETS.items():
        e = ev[sel(ev) & (ev.force_pn >= 18)]
        fi = e[e.category == "intermediate_fwd"]
        nu = un[sel(un)]
        nl = lad[sel(lad)]
        # relative positions vs detector-aware uniform
        for col in ("u", "d1_bp_A", "d1_bp_B", "d2_bp_A", "d2_bp_B", "abs_bp_A", "abs_bp_B"):
            x, y = fi[col].dropna().to_numpy(), nu[col].dropna().to_numpy()
            if x.size < 3 or y.size < 3:
                continue
            ks = ks_2samp(x, y)
            rows.append({"set": sname, "test": "forward intermediate vs detector-aware uniform",
                         "quantity": col, "n_obs": x.size, "n_null": y.size,
                         "median_obs": round(float(np.median(x)), 3), "median_null": round(float(np.median(y)), 3),
                         "stat": round(float(ks.statistic), 3), "p": round(float(ks.pvalue), 4)})
        # landmark proximity (absolute) vs uniform and vs shuffle
        for mdl in ("A", "B"):
            x = fi[f"abs_bp_{mdl}"].dropna().to_numpy()
            if x.size < 3:
                continue
            dn = landmark_dist(nu[f"abs_bp_{mdl}"].dropna().to_numpy())
            dx = landmark_dist(x)
            # shuffle: u permuted across forward intermediates, position recomputed
            sh = []
            u = fi.u.to_numpy()
            l1 = fi[f"prev_long_abs_bp_{mdl}"].to_numpy()
            gp = fi[f"gap_bp_{mdl}"].to_numpy()
            for _ in range(2000):
                sh.append(np.median(landmark_dist(l1 + rng.permutation(u) * gp)))
            sh = np.array(sh)
            rows.append({"set": sname, "test": "distance to nearest 601 landmark (entry, dyad-40, dyad, dyad+40, exit)",
                         "quantity": f"abs_bp_{mdl}", "n_obs": x.size, "n_null": dn.size,
                         "median_obs": round(float(np.median(dx)), 1), "median_null": round(float(np.median(dn)), 1),
                         "stat": round(float(np.median(sh)), 1),
                         "p": round(float(mannwhitneyu(dx, dn, alternative="less").pvalue), 4),
                         "p_shuffle": round(float((np.sum(sh <= np.median(dx)) + 1) / (sh.size + 1)), 4),
                         "note": "stat = median of the u-shuffle null; p = one-sided MWU vs uniform"})
        # comb statistics on d1 (bp) vs detector-aware uniform
        for mdl in ("A", "B"):
            x = fi[f"d1_bp_{mdl}"].dropna().to_numpy()
            y = nu[f"d1_bp_{mdl}"].dropna().to_numpy()
            if x.size < 3 or y.size < 3:
                continue
            for per in (40.0, 80.0):
                cn = np.array([comb_c(rng.choice(y, x.size), per) for _ in range(4000)])
                c_obs = comb_c(x, per)
                h1 = nl[nl["mode"] == f"ladder {per:g} bp {mdl}"][f"d1_bp_{mdl}"].dropna().to_numpy()
                ch = np.array([comb_c(rng.choice(h1, x.size), per) for _ in range(4000)]) if h1.size > 3 else np.array([np.nan])
                rows.append({"set": sname, "test": f"comb C_{per:.0f} = <cos(2 pi d1 / {per:.0f} bp)>",
                             "quantity": f"d1_bp_{mdl}", "n_obs": x.size, "n_null": y.size,
                             "median_obs": round(c_obs, 3), "median_null": round(float(np.median(cn)), 3),
                             "stat": round(float(np.percentile(cn, 95)), 3),
                             "p": round(float((np.sum(cn >= c_obs) + 1) / (cn.size + 1)), 4),
                             "C_ladder_median": round(float(np.nanmedian(ch)), 3),
                             "power_at_5pct": round(float(np.nanmean(ch > np.percentile(cn, 95))), 3),
                             "power_10pct_scale_error": round(float(np.nanmean(np.array(
                                 [comb_c(rng.choice(h1, x.size) * rng.normal(1.0, 0.10, x.size), per)
                                  for _ in range(2000)]) > np.percentile(cn, 95))), 3) if h1.size > 3 else np.nan,
                             "note": "stat = 95th pct of the uniform null C; power = H1 ladder (k x period bp, re-detected)"})
        # excursion amplitudes vs FP excursions (N2, N3) and vs 10 pN real excursions
        ex = e[e.category.str.startswith("excursion") & e.returns]
        ex10 = ev[sel(ev) & (ev.force_pn == 10.0) & ev.category.str.startswith("excursion") & ev.returns]
        for nm_, ref in (("N2 boot FP", fp[(fp.null == "N2 boot") & fp.category.str.startswith("excursion") & fp.returns]),
                         ("N3 creep FP", fp[(fp.null == "N3 creep") & fp.category.str.startswith("excursion") & fp.returns]),
                         ("10 pN real", ex10)):
            ref_amp = (ref.amp_nm if "amp_nm" in ref else ref.amp_prev_nm).abs().to_numpy()
            x = ex.amp_nm.abs().to_numpy()
            if x.size < 3 or ref_amp.size < 3:
                rows.append({"set": sname, "test": f"reversible excursion |amplitude| vs {nm_}",
                             "quantity": "abs amp_nm", "n_obs": x.size, "n_null": ref_amp.size})
                continue
            ks = ks_2samp(x, ref_amp)
            rows.append({"set": sname, "test": f"reversible excursion |amplitude| vs {nm_}",
                         "quantity": "abs amp_nm", "n_obs": x.size, "n_null": ref_amp.size,
                         "median_obs": round(float(np.median(x)), 1), "median_null": round(float(np.median(ref_amp)), 1),
                         "stat": round(float(ks.statistic), 3), "p": round(float(ks.pvalue), 4)})
        for mdl in ("A", "B"):
            x = ex[f"amp_bp_{mdl}"].abs().dropna().to_numpy()
            if x.size < 3:
                continue
            for per in (40.0, 80.0):
                rows.append({"set": sname, "test": f"excursion |amp| comb C_{per:.0f}", "quantity": f"amp_bp_{mdl}",
                             "n_obs": x.size, "median_obs": round(comb_c(x, per), 3),
                             "note": "> 0 = amplitudes near multiples of the period; no detector-aware null"})
    return pd.DataFrame(rows)


def rate_table(ev, pc, cnt, dn, lam, n_boot=2000, seed=5) -> pd.DataFrame:
    """Brief events between long stalls per OPENED plateau (some long stall
    > 100 bp, model A; plateaus where the fork never left the start cannot host
    an event between two main stalls), per force and pooled >= 18 pN, 95 % CI by
    resampling tethers; the same plateaus' matched-null expectations (N1 white,
    N2 boot, N3 creep); forward intermediates supported as discrete levels
    (BIC(ramp) - BIC(steps) > 6) against the same count on N3."""
    rng = np.random.default_rng(seed)
    det = det_name(lam)
    nl = cnt[cnt.detector == det]
    dnp = dn[dn.kind == "plateau"]
    rows = []
    for sname, sel in SETS.items():
        P_ = pc[sel(pc) & (pc.max_long_bp_A > 100)].assign(tether=lambda d: d.data_id.astype(str) + "_" + d.bead)
        E_ = ev[sel(ev) & ev.category.isin(CATS)].assign(tether=lambda d: d.data_id.astype(str) + "_" + d.bead)
        groups = [(f"{f:g}", g) for f, g in P_.groupby("force_pn")] + [(">=18", P_[P_.force_pn >= 18])]
        for f, pg in groups:
            eg = E_[E_.i.isin(pg.i)]
            teth = pg.tether.unique()
            if not teth.size:
                continue
            nper = pg.groupby("tether").size().reindex(teth).to_numpy()
            pick = rng.integers(0, teth.size, (n_boot, teth.size))
            dsub = dnp[dnp.i.isin(pg.i)]
            for cat in ("between (all 4)", "intermediate_fwd", "intermediate_fwd discrete (dBIC>6)",
                        "excursion (returns)"):
                if cat == "between (all 4)":
                    sub, ncol = eg, list(CATS)
                elif cat == "intermediate_fwd":
                    sub, ncol = eg[eg.category == "intermediate_fwd"], ["intermediate_fwd"]
                elif cat.startswith("intermediate_fwd discrete"):
                    sub, ncol = eg[(eg.category == "intermediate_fwd") & (eg.dbic_ramp_minus_steps > DBIC_STRONG)], None
                else:
                    sub, ncol = eg[eg.category.str.startswith("excursion") & eg.returns], ["excursion_returns"]
                c = sub.groupby("tether").size().reindex(teth, fill_value=0).to_numpy()
                boots = c[pick].sum(1) / nper[pick].sum(1)
                exp = {}
                for nn in ("N1 white", "N2 boot", "N3 creep"):
                    if ncol is None:
                        g_ = dsub[dsub.null == nn]
                        exp[nn] = g_.n_discrete.sum() / N_NULL / len(pg) if len(g_) else np.nan
                    else:
                        g_ = nl[(nl.null == nn) & nl.i.isin(pg.i)]
                        exp[nn] = g_[ncol].sum().sum() / N_NULL / len(pg)
                fwd_moves = int(pg.n_fwd_moves.sum())
                rows.append({"set": sname, "force_pn": f, "category": cat, "n_tethers": teth.size,
                             "n_opened_plateaus": len(pg), "n_events": int(c.sum()),
                             "per_opened_plateau": round(c.sum() / len(pg), 4),
                             "ci_lo": round(float(np.percentile(boots, 2.5)), 4),
                             "ci_hi": round(float(np.percentile(boots, 97.5)), 4),
                             "n_major_fwd_moves": fwd_moves,
                             "per_major_fwd_move": round(c.sum() / fwd_moves, 4) if fwd_moves else np.nan,
                             **{f"null_{k.split()[0]}_per_opened_plateau": round(float(v), 4) for k, v in exp.items()}})
    return pd.DataFrame(rows)


def disc_one(args) -> list[dict]:
    """Forward intermediates and their BIC(ramp) - BIC(steps) on the real trace
    and on N2 (noise) / N3 (continuous ramps over the same intervals) nulls."""
    i, lam = args
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    rng = np.random.default_rng(77_000 + i)
    st0 = detect(t, y, sigma_diff(y), LAM_REF, f)
    res = residuals(st0, y)
    sk, cr = skeleton(st0, y.size), creep_base(st0, t, y)
    traces = [("real", 0, y)] + [(nm_, r_, b_ + block_boot(res, y.size, rng))
                                 for r_ in range(N_NULL) for nm_, b_ in (("N2 boot", sk), ("N3 creep", cr))]
    rows = []
    for name, rep, ys in traces:
        s_ = sigma_diff(ys)
        st = detect(t, ys, s_, lam, f)
        lev = st.level_nm.to_numpy()
        longs = st.index[st.long].to_numpy()
        n_int = n_disc = 0
        for i1, i2 in zip(longs[:-1], longs[1:]):
            mid = st.loc[i1 + 1:i2 - 1]
            k = int((mid.category == "intermediate_fwd").sum())
            if not k:
                continue
            db = float(ramp_vs_steps(st, int(i1), int(i2), ys, s_, t))
            n_int += k
            n_disc += k if db > DBIC_STRONG else 0
            rows.append({"i": i, "force_pn": f, "excellent": q["excellent"], "null": name, "rep": rep,
                         "kind": "transition", "n_intermediate_fwd": k, "dbic": db,
                         "gap_nm": lev[i2] - lev[i1]})
        rows.append({"i": i, "force_pn": f, "excellent": q["excellent"], "null": name, "rep": rep,
                     "kind": "plateau", "n_intermediate_fwd": n_int, "n_discrete": n_disc})
    return rows


def stage_disc(plats, qc, lam) -> pd.DataFrame:
    idx = qc[qc.passes].i.tolist()
    res = pool_map(disc_one, [(i, lam) for i in idx], plats, chunksize=4)
    dn = pd.DataFrame([r for rr in res for r in rr])
    dn.to_csv(OUT / "discreteness_ramp_vs_steps.csv.gz", index=False)
    return dn


def disc_summary(dn: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sname, sel in SETS.items():
        d = dn[sel(dn) & (dn.force_pn >= 18)]
        pl = d[d.kind == "plateau"]
        tr_ = d[d.kind == "transition"]
        nplat = pl[pl.null == "real"].i.nunique()
        for nm_ in ("real", "N2 boot", "N3 creep"):
            reps = 1 if nm_ == "real" else N_NULL
            a = pl[pl.null == nm_]
            b = tr_[tr_.null == nm_]
            rows.append({"set": sname, "null": nm_, "n_plateaus": nplat,
                         "fwd_intermediates_per_plateau": round(a.n_intermediate_fwd.sum() / reps / max(nplat, 1), 4),
                         "discrete_dBIC_gt6_per_plateau": round(a.n_discrete.sum() / reps / max(nplat, 1), 4),
                         "transitions_with_intermediate": round(len(b) / reps, 1),
                         "frac_transitions_dBIC_gt6": round(float((b.dbic > DBIC_STRONG).mean()), 3) if len(b) else np.nan,
                         "median_dBIC": round(float(b.dbic.median()), 1) if len(b) else np.nan})
    return pd.DataFrame(rows)


# ================================================================= driver
STAGES = ("prep", "scan", "inject", "events", "disc", "geom", "figs")


def run(stage: str, plats, data) -> int:
    if stage == "all":
        todo = STAGES
    elif stage.endswith("+"):
        todo = STAGES[STAGES.index(stage[:-1]):]
    else:
        todo = (stage,)
    stats = {}
    sp = OUT / "stats.json"
    if sp.exists():
        stats = json.loads(sp.read_text())

    t0 = time.time()
    if "prep" in todo:
        qc = stage_prep(plats)
    else:
        qc = pd.read_csv(OUT / "plateaus_qc.csv")
    g = qc[qc.passes]
    stats["plateaus"] = {"collected": len(qc), "pass": int(len(g)),
                         "fail_reasons": qc[~qc.passes].gate_reason.value_counts().to_dict(),
                         "tethers_pass": int(g.groupby(["data_id", "bead"]).ngroups),
                         "excellent_tethers_pass": int(g[g.excellent].groupby(["data_id", "bead"]).ngroups),
                         "excellent_plateaus_pass": int(g.excellent.sum()),
                         "noise_median_sigma_nm": round(float(g.sigma_nm.median()), 1),
                         "noise_median_acf1_3_5_10": [round(float(g[f"acf{k}"].median()), 3) for k in (1, 3, 5, 10)],
                         "noise_acf1_p95": round(float(g.acf1.quantile(0.95)), 3),
                         "sigma_res_over_sigma_diff": round(float((g.sigma_res_nm / g.sigma_nm).median()), 3),
                         "fps_by_force": g.groupby("force_pn").fps.median().to_dict()}
    print(f"[{time.time() - t0:.0f} s] prep done", flush=True)
    if todo[-1] == "prep":
        sp.write_text(json.dumps(stats, indent=1, default=str))
        return 0

    if "scan" in todo:
        cnt, evn, lam = stage_scan(plats, qc)
        stats["lambda"] = lam
    else:
        cnt = pd.read_csv(OUT / "detector_null_scan_counts.csv.gz")
        evn = pd.read_csv(OUT / "detector_null_scan_events.csv.gz")
        lam = choose_lambda(cnt)
        stats["lambda"] = lam
    summ = scan_summary(cnt)
    stats["detector_scan"] = summ[summ.forces == ">=18 pN"].set_index(["detector", "null"]) \
        .between_per_plateau.unstack().round(4).to_dict()
    stats["detector_scan_10pN"] = summ[summ.forces == "10 pN"].set_index(["detector", "null"]) \
        .between_per_plateau.unstack().round(4).to_dict()
    print(f"[{time.time() - t0:.0f} s] scan done: lambda = {lam:g}", flush=True)
    if todo[-1] == "scan":
        sp.write_text(json.dumps(stats, indent=1, default=str))
        return 0

    if "inject" in todo:
        inj = stage_inject(lam, plats, qc)
    else:
        inj = pd.read_csv(OUT / "detector_injection_recall.csv")
    lim = detection_limits(inj)
    stats["detection_limits_sigma36"] = {
        f"{d} {k} {fps}": {f"{a:g} nm": v for a, v in zip(g.amp_nm, g.shortest_dwell_recall80_s)}
        for (d, k, fps), g in lim[lim.sigma_nm == 36.0].groupby(["detector", "kind", "fps_class"])}
    stats["frame_clock"] = {c: {"mean_frames_in_0.1s": float(inj[(inj.fps_class == c) & (inj.duration_s == 0.1)].mean_frames_in_event.mean()),
                                "frac_0.05s_event_in_gap": float(inj[(inj.fps_class == c) & (inj.duration_s == 0.05)].frac_event_in_gap.mean())}
                            for c in INJ_CLASSES}
    print(f"[{time.time() - t0:.0f} s] inject done", flush=True)
    if todo[-1] == "inject":
        sp.write_text(json.dumps(stats, indent=1, default=str))
        return 0

    if "events" in todo:
        out = stage_events(plats, qc, lam)
        long_threshold_table(out).to_csv(OUT / "long_threshold_sensitivity.csv", index=False)
        ev, lg, tr, pc = out[LONG_S]
    else:
        ev = pd.read_csv(OUT / "events_brief_all.csv")
        lg = pd.read_csv(OUT / "long_stalls.csv")
        tr = pd.read_csv(OUT / "transitions_long_to_long.csv")
        pc = pd.read_csv(OUT / "plateau_counts.csv")
    print(f"[{time.time() - t0:.0f} s] events done", flush=True)
    if todo[-1] == "events":
        sp.write_text(json.dumps(stats, indent=1, default=str))
        return 0

    if "disc" in todo:
        dn = stage_disc(plats, qc, lam)
    else:
        dn = pd.read_csv(OUT / "discreteness_ramp_vs_steps.csv.gz")
    ds_ = disc_summary(dn)
    ds_.to_csv(OUT / "discreteness_summary.csv", index=False)
    stats["discreteness"] = ds_.to_dict(orient="records")
    print(f"[{time.time() - t0:.0f} s] discreteness done", flush=True)
    print(ds_.to_string(index=False), flush=True)
    if todo[-1] == "disc":
        sp.write_text(json.dumps(stats, indent=1, default=str))
        return 0

    if "geom" in todo:
        un = stage_geom(plats, ev, tr, lam)
    else:
        un = pd.read_csv(OUT / "geometry_null_uniform_detector_aware.csv.gz")
    ladder_table(ev[ev.category.isin(CATS)], tr).to_csv(OUT / "transitions_with_intermediates_ladder.csv", index=False)
    gt = geometry_tests(ev, un, evn, lam)
    gt.to_csv(OUT / "geometry_tests.csv", index=False)
    rt = rate_table(ev, pc, cnt, dn, lam)
    rt.to_csv(OUT / "rate_by_force.csv", index=False)
    print(f"[{time.time() - t0:.0f} s] geometry / rates done", flush=True)

    if "figs" in todo or stage == "all":
        make_figures(plats, qc, cnt, evn, inj, ev, lg, tr, pc, un, rt, dn, lam)
        print(f"[{time.time() - t0:.0f} s] figures done", flush=True)
    stats.update(summarise(qc, ev, lg, tr, pc, un, gt, rt, lam))
    sp.write_text(json.dumps(stats, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    return 0


def summarise(qc, ev, lg, tr, pc, un, gt, rt, lam) -> dict:
    s = {}
    for sname, sel in SETS.items():
        e = ev[sel(ev) & (ev.force_pn >= 18)]
        b = e[e.category.isin(CATS)]
        c = pc[sel(pc) & (pc.force_pn >= 18)]
        t_ = tr[sel(tr) & (tr.force_pn >= 18)]
        fi = b[b.category == "intermediate_fwd"]
        ex = b[b.category.str.startswith("excursion") & b.returns]
        s[sname] = {"plateaus_ge18": len(c), "tethers_ge18": int(c.groupby(["data_id", "bead"]).ngroups),
                    "long_stalls_ge18": int(lg[sel(lg) & (lg.force_pn >= 18)].shape[0]),
                    "fwd_long_to_long_moves": int((t_.direction == "fwd").sum()),
                    "fwd_moves_with_intermediate": int(((t_.direction == "fwd") & (t_.n_intermediate > 0)).sum()),
                    "between_events": int(len(b)), "by_category": b.category.value_counts().to_dict(),
                    "excursions_returning": int(len(ex)),
                    "excursions_returning_opened_unflagged": int((~ex.at_closed_fork.astype(bool) & ~ex.flag_tracking
                                                                  & ~ex.flag_large_fast).sum()),
                    "between_events_tracking_flag": int((b.flag_tracking | b.flag_large_fast).sum()),
                    "between_events_at_closed_fork": int(b.at_closed_fork.astype(bool).sum()),
                    "fwd_intermediate_discrete_dBIC_gt6": int((fi.dbic_ramp_minus_steps > DBIC_STRONG).sum()),
                    "fwd_intermediate_median_duration_s": round(float(fi.duration_s.median()), 3) if len(fi) else None,
                    "fwd_intermediate_median_d1_bp_A": round(float(fi.d1_bp_A.median()), 1) if len(fi) else None,
                    "fwd_intermediate_median_d1_bp_B": round(float(fi.d1_bp_B.median()), 1) if len(fi) else None,
                    "fwd_intermediate_median_u": round(float(fi.u.median()), 3) if len(fi) else None,
                    "excursion_median_abs_amp_bp_A": round(float(ex.amp_bp_A.abs().median()), 1) if len(ex) else None,
                    "excursion_median_abs_amp_bp_B": round(float(ex.amp_bp_B.abs().median()), 1) if len(ex) else None,
                    "excursion_median_duration_s": round(float(ex.duration_s.median()), 3) if len(ex) else None}
    return {"summary": s, "lambda": lam}


# ================================================================= figures
CAT_COL = {"intermediate_fwd": "#E64B35", "intermediate_back": "#F39B7F",
           "excursion_fwd": "#3C5488", "excursion_back": "#00A087"}
FPS_COL = {"100 fps": "#E64B35", "50 fps": "#3C5488", "<=33 fps": "#00A087"}


def _export(fig, name):
    from whz_plot import export_panel
    export_panel(fig, OUT / name, formats=("pdf", "png"), dpi=300)


def _shade_sites(ax, horizontal=False):
    span = ax.axhspan if horizontal else ax.axvspan
    line = ax.axhline if horizontal else ax.axvline
    for lo, hi in (N1, N2):
        span(lo, hi, color="#BFBFBF", alpha=0.25, lw=0, zorder=0)
    for dy in DYADS:
        line(dy, color="#666666", lw=0.6, zorder=0)
        for x in (dy - 40, dy + 40):
            line(x, color="#666666", lw=0.5, ls=":", zorder=0)


def fig_validation(qc, cnt, inj):
    from whz_plot import style_axes, fit_line
    g = qc[qc.passes]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.6))
    ax = axes[0, 0]
    lags = np.arange(1, 11)
    acf = g[[f"acf{k}" for k in lags]].to_numpy()
    med = np.median(acf, 0)
    lo, hi = np.percentile(acf, 5, 0), np.percentile(acf, 95, 0)
    ax.errorbar(lags, med, yerr=[med - lo, hi - med], fmt="o", ms=3.5, color="black", lw=0.8, capsize=2)
    ax.axhline(0, color="#BFBFBF", lw=0.8, zorder=0)
    style_axes(ax, xlabel="Lag (frames)", ylabel="Residual autocorrelation", ylim=(-0.3, 0.3))
    ax.set_title("a  frame noise is white (median, 5-95 %)", fontsize=8, loc="left")

    ax = axes[0, 1]
    s = scan_summary(cnt)
    s = s[s.forces == ">=18 pN"]
    cols = {"real": "black", "N1 white": "#8491B4", "N2 boot": "#3C5488", "N3 creep": "#E64B35"}
    for nul, c in cols.items():
        y = [float(s[(s.detector == det_name(l_)) & (s.null == nul)].between_per_plateau.iloc[0]) for l_ in LAMBDAS]
        fit_line(ax, LAMBDAS, np.maximum(y, 1e-4), color=c, linewidth=1.0, label=nul, marker="o", ms=3)
        yb = float(s[(s.detector == det_name("BS")) & (s.null == nul)].between_per_plateau.iloc[0])
        ax.plot([7.2], [max(yb, 1e-4)], marker="s", ms=4, mfc="none", mec=c, ls="none")
    ax.set_yscale("log")
    ax.set_xticks(list(LAMBDAS) + [7.2])
    ax.set_xticklabels([f"{l_:g}" for l_ in LAMBDAS] + ["BS"])
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    style_axes(ax, xlabel="OP penalty lambda", ylabel="Events between long stalls / plateau")
    ax.set_title("b  false events on matched nulls (>= 18 pN)", fontsize=8, loc="left")

    lim = detection_limits(inj)
    for ax, kind, ttl in ((axes[1, 0], "fwd", "c  forward intermediate 0-A-2A"),
                          (axes[1, 1], "exc", "d  reversible excursion 0-A-0")):
        for fps, c in FPS_COL.items():
            for det, ls, lab in ((None, "-", fps), ("BS (dwell_steps)", "--", None)):
                dname = [d for d in lim.detector.unique() if d.startswith("OP")][0] if det is None else det
                q = lim[(lim.fps_class == fps) & (lim.sigma_nm == 36.0) & (lim.detector == dname) & (lim.kind == kind)]
                q = q.dropna(subset=["shortest_dwell_recall80_s"]).sort_values("amp_nm")
                if len(q):
                    fit_line(ax, q.shortest_dwell_recall80_s, q.amp_nm, color=c, linestyle=ls, linewidth=1.0,
                             marker="o", ms=2.5, label=lab)
        ax.plot([], [], color="#666666", ls="--", lw=1.0, label="BS (dwell_steps)")
        ax.set_xscale("log")
        ax.legend(frameon=False, fontsize=7, loc="lower left")
        style_axes(ax, xlabel="Shortest dwell with recall >= 0.8 (s)", ylabel="Amplitude A (nm)",
                   xlim=(0.015, 1.2), ylim=(0, 215))
        ax.set_title(ttl + ", sigma 36 nm", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig1_detector_validation")


def fig_dwell(plats, qc, lam):
    from whz_plot import style_axes
    ev = pd.read_csv(OUT / "events_brief_all.csv")
    lg = pd.read_csv(OUT / "long_stalls.csv")
    sens = pd.read_csv(OUT / "long_threshold_sensitivity.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.8), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    bins = np.logspace(np.log10(0.02), np.log10(16), 34)
    allseg = pd.concat([ev[["force_pn", "duration_s", "censored"]], lg[["force_pn", "duration_s", "censored"]]])
    for sel, c, lab in ((allseg.force_pn >= 18, "black", ">= 18 pN"), (allseg.force_pn == 10.0, "#8491B4", "10 pN")):
        d = allseg[sel & ~allseg.censored].duration_s
        ax.hist(d, bins=bins, histtype="step", color=c, lw=1.0, label=f"{lab}, uncensored")
    for x, ls in ((LONG_S, "-"),) + tuple((x, ":") for x in LONG_SENS):
        ax.axvline(x, color="#666666", ls=ls, lw=0.8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    style_axes(ax, xlabel="Level duration (s)", ylabel="Levels per bin")
    ax.set_title("a  all fitted levels, OP lambda = %g" % lam, fontsize=8, loc="left")
    ax = axes[1]
    for sname, mk in (("all_clean", "o"), ("excellent", "s")):
        q = sens[(sens.set == sname) & (sens.forces == ">=18 pN")].sort_values("long_s")
        ax.plot(q.long_s, q.between_per_plateau, marker=mk, ls="none", color=SET_COL[sname], ms=4,
                label=sname.replace("_", " "))
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 2])
    ax.set_xticklabels(["0.5", "1", "2"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.legend(frameon=False, fontsize=7)
    style_axes(ax, xlabel="Long-stall threshold (s)", ylabel="Brief events between / plateau")
    ax.set_title("b  threshold sensitivity", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig2_dwell_threshold")


def pick_examples(ev, tr, n_int=7, n_exc=3):
    t = tr[(tr.direction == "fwd") & (tr.n_intermediate > 0) & (tr.force_pn >= 18)].copy()
    t["key"] = t.data_id.astype(str) + t.bead
    out, seen, perf = [], set(), {}
    exc_first = t[t.excellent]
    for r in list(exc_first.itertuples()) + list(t.sort_values("dbic_ramp_minus_steps", ascending=False).itertuples()):
        if r.key in seen or perf.get(r.force_pn, 0) >= 2 or r.sigma_nm > 40:
            continue
        out.append(("int", r.i, r.t_l1_end_s, r.t_l2_start_s))
        seen.add(r.key)
        perf[r.force_pn] = perf.get(r.force_pn, 0) + 1
        if len(out) >= n_int:
            break
    ex = ev[ev.category.str.startswith("excursion") & ev.returns & (ev.force_pn >= 18)].copy()
    ex["z"] = ex.amp_nm.abs() / ex.sem_nm
    ex["flag"] = ex.flag_tracking | ex.flag_large_fast | ex.at_closed_fork.astype(bool)
    k = 0
    for r in pd.concat([ex[~ex.flag].sort_values("z", ascending=False), ex[ex.flag].sort_values("z", ascending=False)]).itertuples():
        key = f"{r.data_id}{r.bead}"
        if key in seen:
            continue
        out.append(("exc", r.i, r.prev_long_end_s, r.next_long_start_s))
        seen.add(key)
        k += 1
        if k >= n_exc:
            break
    return out


def fig_examples(plats, ev, tr, lam):
    from whz_plot import style_axes, raw_trace, fit_line
    ex = pick_examples(ev, tr)
    n = len(ex)
    nrow = int(np.ceil(n / 2))
    fig, axes = plt.subplots(nrow, 2, figsize=(7.09, 1.85 * nrow + 0.3), squeeze=False)
    rows = []
    for ax, (kind, i, ta, tb) in zip(axes.flat, ex):
        q = plats[i]
        y, t, f = q["y"], q["t"], q["force_pn"]
        sig = sigma_diff(y)
        st = detect(t, y, sig, lam, f)
        ga, gb = gains(f)
        lo, hi = max(0.0, ta - 1.5), min(t[-1], tb + 1.5)
        w = (t >= lo) & (t <= hi)
        raw_trace(ax, t[w], y[w], color="#666666", filled=True, marker_area=2.0, edge_width=0.0, alpha=0.7,
                  rasterized=True)
        xs, ys = [], []
        for r in st.itertuples():
            a_, b_ = max(r.t0, lo), min(r.t1, hi)
            if b_ > a_:
                xs += [a_, b_]
                ys += [r.level_nm, r.level_nm]
        fit_line(ax, xs, ys, color="black", linewidth=0.9)
        sub = ev[(ev.i == i) & ev.category.isin(CATS) & (ev.t_start_s >= ta - 1e-3) & (ev.t_start_s < tb)]
        for r in sub.itertuples():
            fit_line(ax, [r.t_start_s, r.t_start_s + r.duration_s], [r.level_nm] * 2, color=CAT_COL[r.category],
                     linewidth=2.2, zorder=5)
        title = f"D{q['data_id']:03d} {q['bead']} step {q['step']} {f:g} pN"
        if kind == "int":
            i1 = int(np.flatnonzero(st.long.to_numpy() & np.isclose(st.t1.to_numpy(), ta, atol=2e-3))[0])
            i2 = int(np.flatnonzero(st.long.to_numpy() & np.isclose(st.t0.to_numpy(), tb, atol=2e-3))[0])
            db, (a_, b_, rf) = ramp_vs_steps(st, i1, i2, y, sig, t, return_fit=True)
            fit_line(ax, t[a_:b_], rf, color="#4DBBD5", linewidth=0.9, linestyle="--", zorder=4)
            l1 = st.level_nm[i1]
            for bp in (40, 80, 120, 160):
                ax.axhline(l1 + bp * ga, color="#BFBFBF", ls=":", lw=0.6, zorder=0)
            title += f", dBIC {db:+.0f}"
        else:
            r0 = sub.iloc[0] if len(sub) else None
            if r0 is not None:
                title += f", {r0.amp_bp_A:+.0f} bp A / {r0.amp_bp_B:+.0f} bp B"
                if r0.flag_tracking or r0.flag_large_fast:
                    title += ", tracking flag"
        for r in sub.itertuples():
            rows.append({"panel": len(rows), "kind": kind, "data_id": q["data_id"], "bead": q["bead"],
                         "step": q["step"], "force_pn": f, "category": r.category, "t_start_s": r.t_start_s,
                         "duration_s": r.duration_s, "d1_bp_A": getattr(r, "d1_bp_A"), "d1_bp_B": r.d1_bp_B,
                         "d2_bp_A": r.d2_bp_A, "d2_bp_B": r.d2_bp_B, "abs_bp_A": r.abs_bp_A, "abs_bp_B": r.abs_bp_B,
                         "dbic": r.dbic_ramp_minus_steps})
        ax.set_xlim(lo, hi)
        ax.set_title(title, fontsize=6.5, loc="left", pad=2)
        style_axes(ax, xlabel="Time in plateau (s)", ylabel="Extension (nm)")
        ax.tick_params(labelsize=6.5)
        ax.xaxis.label.set_size(7)
        ax.yaxis.label.set_size(7)
    for ax in list(axes.flat)[n:]:
        ax.set_visible(False)
        ax.remove()
    fig.tight_layout(h_pad=0.6, w_pad=0.8)
    _export(fig, "fig3_examples")
    pd.DataFrame(rows).to_csv(OUT / "fig3_examples_events.csv", index=False)


def fig_relative(ev, un, evn, lam):
    from whz_plot import style_axes
    det = det_name(lam)
    nu = un[un["mode"] == "uniform"]
    e = ev[(ev.force_pn >= 18)]
    fi = e[e.category == "intermediate_fwd"]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.2))

    def panel(ax, x, xn, bins, xl, ttl, guides=()):
        ax.hist(x, bins=bins, color="#E64B35", alpha=0.75, lw=0, label=f"observed (n = {len(x)})")
        xe = fi[fi.excellent][xl[1]].dropna() if xl[1] in fi else pd.Series(dtype=float)
        if len(xe):
            ax.hist(xe, bins=bins, color="black", histtype="step", lw=1.0, label=f"excellent (n = {len(xe)})")
        if len(xn):
            h, b = np.histogram(xn, bins=bins)
            ax.stairs(h * len(x) / max(h.sum(), 1), b, color="black", ls="--", lw=1.0, label="detector-aware uniform")
        for gx in guides:
            ax.axvline(gx, color="#BFBFBF", ls=":", lw=0.7, zorder=0)
        ax.legend(frameon=False, fontsize=6.5)
        style_axes(ax, xlabel=xl[0], ylabel="Forward intermediates")
        ax.set_title(ttl, fontsize=8, loc="left")

    panel(axes[0, 0], fi.u.dropna(), nu.u.dropna(), np.linspace(0, 1, 11),
          ("Fraction of the gap u = d1 / (L2 - L1)", "u"), "a  position between the two long stalls")
    panel(axes[0, 1], fi.d1_bp_A.dropna(), nu.d1_bp_A.dropna(), np.arange(0, 720, 20),
          ("d1 = intermediate - preceding stall (bp, model A)", "d1_bp_A"), "b  d1, model A (2xss - xds)",
          guides=np.arange(40, 720, 40))
    panel(axes[1, 0], fi.d1_bp_B.dropna(), nu.d1_bp_B.dropna(), np.arange(0, 500, 20),
          ("d1 = intermediate - preceding stall (bp, model B)", "d1_bp_B"), "c  d1, model B (2xss)",
          guides=np.arange(40, 500, 40))
    ax = axes[1, 1]
    ex = e[e.category.str.startswith("excursion") & e.returns]
    bins = np.arange(0, 620, 25)
    ok = ex[~(ex.flag_tracking | ex.flag_large_fast)]
    ax.hist(ok.amp_bp_A.abs(), bins=bins, color="#3C5488", alpha=0.8, lw=0, label=f"observed, unflagged (n = {len(ok)})")
    ax.hist(ex.amp_bp_A.abs(), bins=bins, color="#3C5488", histtype="step", lw=1.0, label=f"observed, all (n = {len(ex)})")
    fp = evn[(evn.detector == det) & (evn.force_pn >= 18) & evn.category.str.startswith("excursion") & evn.returns]
    for nn, ls in (("N2 boot", "--"), ("N3 creep", ":")):
        a = fp[fp.null == nn]
        amp = (a.amp_prev_nm.abs() / g_a(a.force_pn.to_numpy()))
        h, b = np.histogram(amp, bins=bins)
        ax.stairs(h / N_NULL, b, color="black", ls=ls, lw=1.0, label=f"{nn} false (per data set)")
    ax.legend(frameon=False, fontsize=6.5)
    style_axes(ax, xlabel="|Excursion amplitude| (bp, model A)", ylabel="Reversible excursions")
    ax.set_title("d  reversible excursions vs matched nulls", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig4_relative_positions")


def fig_absolute(ev, lg, un):
    from whz_plot import style_axes
    nu = un[un["mode"] == "uniform"]
    e = ev[(ev.force_pn >= 18)]
    fi = e[e.category == "intermediate_fwd"]
    ex = e[e.category.str.startswith("excursion") & e.returns]
    L = lg[(lg.force_pn >= 18)]
    fig, axes = plt.subplots(2, 1, figsize=(7.09, 4.6))
    for ax, mdl, top in ((axes[0], "A", 1200), (axes[1], "B", 850)):
        bins = np.arange(-100, top + 1, 20)
        _shade_sites(ax)
        h, b = np.histogram(L[f"abs_bp_{mdl}"], bins=bins)
        ax.stairs(h / h.sum() / 20, b, color="#666666", lw=0.9, label=f"long stalls (n = {len(L)})")
        h, b = np.histogram(fi[f"abs_bp_{mdl}"], bins=bins)
        ax.stairs(h / max(h.sum(), 1) / 20, b, color="#E64B35", fill=True, alpha=0.6,
                  label=f"forward intermediates (n = {len(fi)})")
        h, b = np.histogram(nu[f"abs_bp_{mdl}"], bins=bins)
        ax.stairs(h / max(h.sum(), 1) / 20, b, color="black", ls="--", lw=0.9, label="uniform null")
        ax.plot(ex[f"abs_bp_{mdl}"], np.full(len(ex), -0.0004), "|", color="#3C5488", ms=6,
                label=f"reversible excursions (n = {len(ex)})")
        ax.axvline(980, color="black", lw=0.6, ls="-.")
        ax.legend(frameon=False, fontsize=6.5, loc="upper right")
        style_axes(ax, xlabel=f"Absolute Luz, model {mdl} (bp)", ylabel="Density (1/bp)", xlim=(-100, top))
    axes[0].set_title("a  model A (2xss - xds); grey: 601 sites, solid dyad, dotted dyad +-40", fontsize=8, loc="left")
    axes[1].set_title("b  model B (2xss); same map", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig5_absolute_positions")


def fig_rates(rt):
    from whz_plot import style_axes
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.9))
    r = rt[rt.force_pn != ">=18"].copy()
    r["F"] = r.force_pn.astype(float)
    for ax, cat, ttl in ((axes[0], "between (all 4)", "a  all brief events between long stalls"),
                         (axes[1], "intermediate_fwd discrete (dBIC>6)", "b  discrete forward intermediates")):
        for sname, mk, dx in (("all_clean", "o", 0.0), ("excellent", "s", 0.35)):
            q = r[(r.set == sname) & (r.category == cat)].sort_values("F")
            ax.errorbar(q.F + dx, q.per_opened_plateau, yerr=[q.per_opened_plateau - q.ci_lo, q.ci_hi - q.per_opened_plateau],
                        fmt=mk, ms=3.5, color=SET_COL[sname], lw=0.8, capsize=2, label=sname.replace("_", " "))
        q = r[(r.set == "all_clean") & (r.category == cat)].sort_values("F")
        ax.plot(q.F, q.null_N2_per_opened_plateau, color="black", ls="--", lw=0.9, label="N2 noise null")
        ax.plot(q.F, q.null_N3_per_opened_plateau, color="black", ls=":", lw=0.9, label="N3 ramp null")
        ax.legend(frameon=False, fontsize=6.5, loc="upper left")
        style_axes(ax, xlabel="Plateau force (pN)", ylabel="Events / opened plateau")
        ax.set_title(ttl, fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig6_rate_vs_force")


def make_figures(plats, qc, cnt, evn, inj, ev, lg, tr, pc, un, rt, dn, lam):
    from whz_plot import apply_style
    apply_style()
    for fn, a in ((fig_validation, (qc, cnt, inj)), (fig_dwell, (plats, qc, lam)),
                  (fig_examples, (plats, ev, tr, lam)), (fig_relative, (ev, un, evn, lam)),
                  (fig_absolute, (ev, lg, un)), (fig_rates, (rt,))):
        fn(*a)
        print(f"  {fn.__name__} done", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recollect", action="store_true")
    ap.add_argument("--stage", default="all", help="all | <stage> | <stage>+ (that stage and the rest)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data = load_cache(args.recollect)
    plats = data["plateaus"]
    _init(plats)
    print(f"{len(plats)} plateaus loaded ({time.time() - t0:.0f} s)", flush=True)
    return run(args.stage, plats, data)


if __name__ == "__main__":
    raise SystemExit(main())
