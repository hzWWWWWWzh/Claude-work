#!/usr/bin/env python3
"""20261009: small jumps (< 100 bp, Luz model B) next to big jumps, per force.

User (2026-10-11): "我希望你能分别统计每一个力的luz peakmap，再zoom进去大跳跃附近的
小跳跃（小于100bp）观察这些小跳变出现的概率" (screenshot: data_026 R039, 10 s
plateaus that sit flat and step up by a few tens of nm near their end). The
per-force Luz peak maps are made elsewhere; this script does the small jumps.

Batch 20261009 only (YTrunk_2601_980bp, 971 bp trunk; 601 sites 314-460, dyad
387, and 511-657, dyad 584). Each titration file holds the force ladder; every
plateau starts from a fork re-zipped during the 1.6 pN rest.

Upstream, read only (nothing upstream is edited or overwritten):
  process_20261009                 load, segment_steps, titration_files, RAW, DEST
  ladder_norm1028_20261009         bead_scale (zero = max extension of the 10 pN
                                   step, rest-dwell drift, usable_steps, t_loss),
                                   MODEL = B (2 x_ss per opened bp)
  stall_structure_briefjumps_20261009 (validated detector, imported as is):
                                   op_segment (exact optimal partitioning,
                                   lambda = 4), sigma_diff, spike_mask,
                                   block_boot, seg_table, handle_nm, g_b
  tracking_quality                 file_tracking

Signal: raw Z at the full frame rate minus the rest-dwell drift line minus the
10 pN zero (nm); single-frame spikes (> 6 robust sigma off a +-3 frame running
median) dropped (same as the validated detector). Luz_B = (ext - zero - 1800 bp
handle stretch from 10 pN) / (2 x_ss(F)).

Plateau gate (same criteria as the validated detector, but the LUT guard in
model B): frame noise 5-60 nm, LUT flicker <= 1 %, clipped <= 2 %, every level
of the lambda = 3 reference fit inside -150..1150 bp B. (The validated run used
-150..1150 bp of model A; a fully unzipped 971 bp fork is ~1400 bp A, so that
guard removed fully opened plateaus; model B fixes that.) Only steps in the
bead's usable_steps (survived; nothing after a loss / fly-away).

Settling: the first T_SKIP = 0.3 s of each plateau are not searched (justified
in stage 'settle': change-point density and the 10 pN control vs time since
arrival). Arrival (end of the 2 mm/s up move, t = 0) is a reference big event.

Events (OP lambda = 4 on t >= T_SKIP): every change point is a step;
  BIG    forward step >= 100 bp B;
  SMALL  |step| < 100 bp B (forward / backward separately);
  big backward (<= -100 bp B) listed, not used as a reference.
Open-fork time: time on levels with 50 <= Luz_B < 900 bp (fork left the start,
not fully unzipped). Closed (< 50 bp) and fully open (>= 900 bp) levels and the
10 pN plateaus are internal no-fork controls.

Nulls (per plateau, its own clock, 3 replicas):
  N2  flat big-jump skeleton (OP levels merged across small steps) + circular
      block bootstrap (10 frames) of the residuals: noise false positives;
  N3  each stretch between big events replaced by a CONTINUOUS piecewise-linear
      curve with 0-3 knots (BIC; knots on a >= 0.1 s grid, i.e. ramps >= 0.1 s)
      fitted to the real data, big jumps kept, + bootstrap residuals: creep cut
      into pieces (reviewer's recommendation, replaces the single-ramp null).
Discreteness per small jump: dBIC = BIC(best continuous PL, 0-3 free knots) -
BIC(two flat levels + one jump) in a window of the two flanking levels (+-1 s);
calibrated on N3 detections (continuous truth) and on injected abrupt steps (H1).
Injection (stage inject): single abrupt steps of 10-120 bp B inserted at random
times into each plateau's own N2 base -> recall per plateau, force, amplitude,
flanking dwell and distance to the nearest big event (detection limits A80/A50).
Real-noise calibration (stage injreal): abrupt steps and 0.3 / 1.0 s ramps added
to the REAL traces (real drift, low-frequency noise, own fork events) -> power of
the discreteness test and its false 'abrupt' rate for creep; abrupt 20-80 nm
steps in the 10 pN plateaus -> recall for the no-fork control.
Empirical no-fork controls: 10 pN plateaus (closed fork, below the ~15 pN
unzipping force) and >= 18 pN levels >= 900 bp B, compared amplitude-matched
(nm) and recall-corrected. Copy check: overlapping ROIs (< 60 px, both beads
present): a small jump at the same frame (+-1) as any change point of the
neighbour is a copy suspect and dropped; same-frame coincidences with all other
beads of the file are counted against a time-shift null (sensitivity table).

    .venv/bin/python smalljumps_modelB_20261009.py [--recollect] [--stage all|qc|settle|detect|nulls|inject|injreal|stats|figs|<stage>+]
Outputs (Proecessed Traces/WHZ_20261009/smalljumps_modelB/):
  plateaus_qc.csv, overlapping_roi_pairs.csv, settle_changepoint_density.csv
  steps_all(_flagged).csv.gz       every change point (class, sizes, levels, dBIC, flags)
  small_jumps_open_fork.csv        the small jumps used (>= 18 pN, 50-900 bp B, no copy suspects)
  small_jump_probability_by_force.csv   P(plateau >= 1), rates, N2/N3, detectable range
  detection_recall_by_force.csv, detection_limits_by_force.csv
  controls_no_fork.csv, controls_amplitude_matched.csv
  per_big_event(_windows).csv, enrichment_by_window.csv, rate_vs_dt_to_big.csv
  discreteness_summary.csv, discreteness_calibration.csv, discreteness_*events*
  amplitude_tests.csv, periodicity_scan.csv, position_tests.csv,
  position_occupancy_null.csv, relative_to_landing_level.csv,
  sensitivity_coincidence_overlap.csv, null_*.csv.gz, injection_*.csv.gz, stats.json
  fig1_examples (+ _events.csv), fig2_settling_and_detection, fig3_probability_vs_force,
  fig4_amplitudes, fig5_timing_vs_big_jumps, fig6_positions_601_map, fig6b_relative_to_landing,
  fig7_discreteness
"""
from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import os
import pickle
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.ndimage import median_filter  # noqa: E402

os.environ.setdefault("MODEL", "B")
import process_20261009 as p  # noqa: E402  (sets Agg)
import matplotlib.pyplot as plt  # noqa: E402
import ladder_norm1028_20261009 as LN  # noqa: E402
import stall_structure_briefjumps_20261009 as SB  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

assert LN.MODEL == "B", "Luz model B required (user's choice 2026-10-11)"

OUT = p.DEST / "smalljumps_modelB"
CACHE = OUT / "_cache_plateaus.pkl"
SEL_CSV = p.DEST / "selection_96_tethers.csv"

LAM = 4.0                 # validated OP penalty (noise-null FP <= 10 % of real counts)
LAM_REF = 3.0             # reference fit for the gate / residuals (as validated)
MMIN = SB.MMIN
T_SKIP = 0.3              # s after arrival not searched (stage 'settle')
BIG_BP = 100.0
OPEN_LO, OPEN_HI = 50.0, 900.0
WIN_S = 1.0               # 'next to a big jump'
DISC_WIN_S = 3.0          # discreteness test: the two flanking levels, clipped to +-3 s
FAR_S = 2.0               # baseline: > 2 s from every big event
SIG_MIN_NM, SIG_MAX_NM = 5.0, 60.0
FLICKER_MAX, REPEAT_MAX = 0.01, 0.02
GUARD_BP_B = (-150.0, 1150.0)
OVERLAP_PX = 60.0
N_NULL = 3
KNOT_GRID_S = 0.1         # N3: ramps >= 0.1 s
PL_MAX_KNOTS = 3
PL_MAX_CAND = 40          # per-event test: candidate knots (frame-gap midpoints)
INJ_AMPS_BP = (10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0)
INJ_REPS = 3
PROCS = 9

N1, N2 = (314.0, 460.0), (511.0, 657.0)
DYADS = (387.0, 584.0)
LANDMARKS = SB.LANDMARKS
LM = np.array(sorted(LANDMARKS.values()), float)
SETS = ("user96", "other_clean")
SET_LABEL = {"user96": "your 96", "other_clean": "other clean"}
SET_COL = {"user96": "#E64B35", "other_clean": "#3C5488"}


# ----------------------------------------------------------------- model B
def gain_b(f) -> float:
    return float(SB.g_b(f)[0])


def luz_b(nm, f) -> np.ndarray:
    return SB.abs_bp(nm, f, "B")


# ----------------------------------------------------------------- collect
def user_set() -> set[tuple[int, str]]:
    s = pd.read_csv(SEL_CSV)
    return set(zip(s.data_id.astype(int), s.bead.astype(str)))


def collect() -> dict:
    sel = user_set()
    titr, _ = p.titration_files()
    plats, beads, pairs = [], [], []
    for d in titr:
        run = p.load(d)
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{d:03d}.h5")
        t_abs = run["time_s"]
        t = t_abs - t_abs[0]
        xy = run["xy"]
        present = [r for r in run["rois"] if trk.get(r, {}).get("status") in ("tracked", "lost")]
        for a, b in itertools.combinations(run["rois"], 2):
            if a in xy and b in xy:
                dist = float(np.hypot(xy[a][0] - xy[b][0], xy[a][1] - xy[b][1]))
                if dist < OVERLAP_PX:
                    pairs.append({"data_id": d, "roi_a": f"R{a}", "roi_b": f"R{b}", "dist_px": round(dist, 1),
                                  "a_present": a in present, "b_present": b in present})
        for roi in run["rois"]:
            bead = f"R{roi}"
            q = trk.get(roi, {})
            rec = {"data_id": d, "bead": bead, "tracking": q.get("status", "unknown"),
                   "user96": (d, bead) in sel, "level_jump_um": q.get("level_jump_um", np.nan)}
            beads.append(rec)
            # the user's hand-picked tethers are kept whatever the file-wide tracking
            # verdict (8 of 96 carry a 'level_jump_artefact' somewhere in the file);
            # their bad plateaus are removed by the per-plateau gate (LUT guard)
            ok_status = rec["tracking"] in ("tracked", "lost") or (rec["user96"] and
                                                                  rec["tracking"] == "level_jump_artefact")
            if not ok_status or not steps:
                rec["reason"] = f"tracking {rec['tracking']}"
                continue
            if abs(steps[0]["force_pn"] - 10.0) > 0.5:
                rec["reason"] = "first step not 10 pN"
                continue
            z_raw = run["z_um"][roi] * 1000.0
            z_s = SB.F.smoothed_nm(run["z_um"][roi])
            ext_s, sc = LN.bead_scale(run, steps, roi, q)       # sets run['z_um_raw'] itself
            zero = sc["zero_nm"]
            rec.update(zero_nm=zero, usable_steps=" ".join(map(str, sc["usable_steps"])),
                       skipped=sc["skipped"], drift=sc["drift"])
            if not np.isfinite(zero):
                rec["reason"] = "no 10 pN zero"
                continue
            rec["reason"] = ""
            ext_raw = z_raw - (z_s - ext_s)                    # raw Z minus the rest-dwell drift line
            t_loss = sc["t_loss"] - t_abs[0] if np.isfinite(sc["t_loss"]) else np.inf
            for st in steps:
                pl = st["plateau"]
                if st["step"] not in sc["usable_steps"] or t[pl[-1]] >= t_loss:
                    continue
                y = ext_raw[pl] - zero
                tt = t[pl]
                ok = np.isfinite(y)
                if ok.sum() < 50:
                    continue
                y, tt = y[ok], tt[ok]
                keep = SB.spike_mask(y)
                plats.append({"data_id": d, "bead": bead, "step": st["step"], "force_pn": round(st["force_pn"], 1),
                              "user96": (d, bead) in sel, "tracking": rec["tracking"],
                              "t": (tt[keep] - t[pl[0]]).astype(float), "y": y[keep].astype(float),
                              "t_file_s": float(t[pl[0]]), "n_spikes": int((~keep).sum()),
                              "hold_s": float(tt[-1] - t[pl[0]]),
                              "fps_median": float(1.0 / np.median(np.diff(tt))),
                              "fps_mean": float((tt.size - 1) / (tt[-1] - tt[0])),
                              "max_gap_s": float(np.max(np.diff(tt))),
                              "t_up_s": float(t[pl[0]] - t[st["up"][0]]) if st["up"].size else np.nan,
                              "zero_nm": zero})
        print(f"data_{d:03d}: {len(run['rois'])} ROIs, {sum(1 for q in plats if q['data_id'] == d)} plateaus",
              flush=True)
    return {"plateaus": plats, "beads": pd.DataFrame(beads), "pairs": pd.DataFrame(pairs)}


def load_cache(recollect=False) -> dict:
    if CACHE.exists() and not recollect:
        with open(CACHE, "rb") as fh:
            return pickle.load(fh)
    OUT.mkdir(parents=True, exist_ok=True)
    data = collect()
    with open(CACHE, "wb") as fh:
        pickle.dump(data, fh)
    return data


# ----------------------------------------------------------------- pool
_P: list = []
_X: dict = {}


def _init(plats, extra=None):
    global _P, _X
    _P = plats
    _X = extra or {}


def pool_map(fn, items, plats, extra=None, chunksize=4):
    ctx = mp.get_context("fork")
    with ctx.Pool(PROCS, initializer=_init, initargs=(plats, extra)) as pool:
        return pool.map(fn, items, chunksize=chunksize)


# ----------------------------------------------------------------- fits
def fit_levels(t, y, sig, lam) -> pd.DataFrame:
    return SB.seg_table(t, y, SB.op_segment(y, sig, lam), sig)


def fit_curve(st: pd.DataFrame, n: int) -> np.ndarray:
    out = np.empty(n)
    for r in st.itertuples():
        out[r.a:r.b] = r.level_nm
    return out


def search(q) -> tuple[np.ndarray, np.ndarray]:
    m = q["t"] >= T_SKIP
    return q["t"][m], q["y"][m]


def step_table(st: pd.DataFrame, t: np.ndarray, f: float) -> pd.DataFrame:
    """One row per change point of the OP fit: size (nm, bp B), the two flanking
    levels (Luz_B), their dwells, class."""
    g = gain_b(f)
    lev, sem = st.level_nm.to_numpy(), st.sem_nm.to_numpy()
    rows = []
    for k in range(len(st) - 1):
        b_prev, a_next = int(st.b[k]), int(st.a[k + 1])
        d_nm = lev[k + 1] - lev[k]
        d_bp = d_nm / g
        cls = ("big_fwd" if d_bp >= BIG_BP else "big_back" if d_bp <= -BIG_BP
               else "small_fwd" if d_bp > 0 else "small_back")
        rows.append({"k": k, "idx": a_next, "t_s": 0.5 * (t[b_prev - 1] + t[a_next]),
                     "gap_s": t[a_next] - t[b_prev - 1], "step_nm": d_nm, "step_bp": d_bp,
                     "se_bp": float(np.hypot(sem[k], sem[k + 1]) / g),
                     "from_nm": lev[k], "to_nm": lev[k + 1],
                     "from_bp": float(luz_b(lev[k], f)[0]), "to_bp": float(luz_b(lev[k + 1], f)[0]),
                     "dwell_before_s": float(st.duration_s[k]), "dwell_after_s": float(st.duration_s[k + 1]),
                     "n_before": int(st.n_frames[k]), "n_after": int(st.n_frames[k + 1]),
                     "before_censored": bool(st.cens_start[k]), "after_censored": bool(st.cens_end[k + 1]),
                     "cls": cls})
    return pd.DataFrame(rows, columns=["k", "idx", "t_s", "gap_s", "step_nm", "step_bp", "se_bp", "from_nm",
                                       "to_nm", "from_bp", "to_bp", "dwell_before_s", "dwell_after_s",
                                       "n_before", "n_after", "before_censored", "after_censored", "cls"])


def big_events(stp: pd.DataFrame) -> np.ndarray:
    """Reference big events: arrival (t = 0) and every forward step >= 100 bp B."""
    return np.concatenate([[0.0], stp.loc[stp.cls == "big_fwd", "t_s"].to_numpy(float)])


def relate(stp: pd.DataFrame, big_t: np.ndarray, st: pd.DataFrame, f: float) -> pd.DataFrame:
    """Signed time to the nearest big event, time since the previous / to the next
    one, and positions relative to the level the previous big event landed on."""
    stp = stp.copy()
    lev = st.level_nm.to_numpy()
    big_k = stp.loc[stp.cls == "big_fwd", "k"].to_numpy(int)
    land_lv = np.concatenate([[lev[0]], lev[big_k + 1]])          # level after arrival / after each big jump
    g = gain_b(f)
    dt_near, dt_prev, dt_next, land, near_kind = [], [], [], [], []
    for r in stp.itertuples():
        others = big_t if r.cls != "big_fwd" else big_t[~np.isclose(big_t, r.t_s)]
        d = r.t_s - others
        j = int(np.argmin(np.abs(d)))
        dt_near.append(float(d[j]))
        near_kind.append("arrival" if others[j] == 0.0 else "big jump")
        prev = d[d > 0]
        nxt = d[d < 0]
        dt_prev.append(float(prev.min()) if prev.size else np.nan)
        dt_next.append(float(-nxt.max()) if nxt.size else np.nan)
        jprev = np.flatnonzero(big_t < r.t_s - 1e-9)
        land.append(land_lv[jprev[-1]] if jprev.size else np.nan)
    stp["dt_nearest_big_s"] = dt_near
    stp["nearest_big_kind"] = near_kind
    stp["dt_since_prev_big_s"] = dt_prev
    stp["dt_to_next_big_s"] = dt_next
    land = np.array(land, float)
    stp["landed_bp"] = luz_b(land, f) if len(stp) else []
    stp["from_rel_landed_bp"] = (stp.from_nm - land) / g
    stp["to_rel_landed_bp"] = (stp.to_nm - land) / g
    return stp


# ----------------------------------------------------------------- continuous piecewise-linear fits
_COMBO_CACHE: dict = {}


def _combos(m: int, k: int) -> np.ndarray:
    key = (m, k)
    if key not in _COMBO_CACHE:
        _COMBO_CACHE[key] = (np.zeros((1, 0), int) if k == 0 else
                             np.array(list(itertools.combinations(range(m), k)), int).reshape(-1, k))
    return _COMBO_CACHE[key]


def pl_fit(t, y, knots_cand, sig, max_k=PL_MAX_KNOTS, want_curve=False):
    """Best CONTINUOUS piecewise-linear fit with 0..max_k knots taken from
    knots_cand (exhaustive), BIC = SSE / sig^2 + (2 + 2k) ln n. Returns
    (bic, k, sse, knots, curve)."""
    n = y.size
    tm = float(t.mean())
    tc = t - tm
    yc = y - y.mean()
    kc = np.asarray(knots_cand, float) - tm
    X = np.column_stack([np.ones(n), tc] + [np.maximum(tc - g_, 0.0) for g_ in kc])
    with np.errstate(all="ignore"):      # spurious matmul warnings (numpy 2.0 + Accelerate)
        M = X.T @ X
        b = X.T @ yc
    yy = float(yc @ yc)
    if not (np.isfinite(M).all() and np.isfinite(b).all()):
        M = np.einsum("ni,nj->ij", X, X)
        b = np.einsum("ni,n->i", X, yc)
    best = (np.inf, 0, yy, np.array([]), None)
    ln = np.log(n)
    for k in range(0, max_k + 1):
        if k > kc.size or 2 + 2 * k >= n:
            break
        cb = _combos(kc.size, k)
        idx = np.concatenate([np.tile([0, 1], (len(cb), 1)), cb + 2], axis=1)
        Ms = M[idx[:, :, None], idx[:, None, :]]
        bs = b[idx]
        ridge = 1e-9 * (np.trace(Ms, axis1=1, axis2=2)[:, None, None] + 1.0) * np.eye(k + 2)[None]
        try:
            coef = np.linalg.solve(Ms + ridge, bs[..., None])[..., 0]
        except np.linalg.LinAlgError:
            continue
        sse = yy - np.einsum("ij,ij->i", bs, coef)
        j = int(np.argmin(sse))
        bic = max(float(sse[j]), 0.0) / sig ** 2 + (2 + 2 * k) * ln
        if bic < best[0]:
            curve = None
            if want_curve:
                curve = np.einsum("ni,i->n", X[:, idx[j]], coef[j]) + y.mean()
            best = (bic, k, float(sse[j]), kc[cb[j]] + tm, curve)
    return best


def step_bic(y, sig, mmin=MMIN):
    """Two flat levels with one jump at the best split (3 parameters)."""
    n = y.size
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y * y)])
    s = np.arange(mmin, n - mmin + 1)
    if s.size == 0:
        return np.inf, -1
    sse1 = cs2[s] - cs[s] ** 2 / s
    sse2 = (cs2[n] - cs2[s]) - (cs[n] - cs[s]) ** 2 / (n - s)
    sse = sse1 + sse2
    j = int(np.argmin(sse))
    return float(sse[j]) / sig ** 2 + 3 * np.log(n), int(s[j])


def event_dbic(t, y, st, k, sig, win=DISC_WIN_S):
    """dBIC = BIC(best continuous PL, 0-3 knots at frame times) - BIC(flat step)
    on the window of the two levels flanking change point k (clipped to +-win s;
    a +-1 s window was too short: a single line then beats a clear 25 nm step
    between two multi-second flat levels, D026 R039 25.9 pN).
    > 0: the abrupt step beats every continuous curve with <= 3 knots."""
    a0, b1 = int(st.a[k]), int(st.b[k + 1])
    s_idx = int(st.a[k + 1])
    tc = 0.5 * (t[s_idx - 1] + t[s_idx])
    m = np.arange(a0, b1)
    m = m[(t[m] >= tc - win) & (t[m] <= tc + win)]
    if m.size < 2 * MMIN + 2:
        return np.nan, np.nan, 0
    tw, yw = t[m], y[m]
    cand_idx = np.arange(1, m.size - 1)
    if cand_idx.size > PL_MAX_CAND:
        sloc = int(np.searchsorted(m, s_idx))
        near = np.arange(sloc - 2, sloc + 2)
        near = near[(near >= 1) & (near < m.size - 1)]
        even = np.unique(np.round(np.linspace(1, m.size - 2, PL_MAX_CAND)).astype(int))
        cand_idx = np.unique(np.concatenate([even, near]))
    bic_pl, kpl, _, knots, _ = pl_fit(tw, yw, tw[cand_idx], sig)
    bic_st, _ = step_bic(yw, sig)
    return float(bic_pl - bic_st), int(kpl), int(m.size)


# ----------------------------------------------------------------- stage qc
def qc_one(i: int) -> dict:
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    sig = SB.sigma_diff(y)
    r9 = y - median_filter(y, 9, mode="nearest")
    _, cnt = np.unique(np.round(y, 3), return_counts=True)
    st = fit_levels(t, y, sig, LAM_REF)
    lv = luz_b(st.level_nm.to_numpy(), f)
    rec = {"i": i, "data_id": q["data_id"], "bead": q["bead"], "step": q["step"], "force_pn": f,
           "set": "user96" if q["user96"] else "other_clean", "tracking": q["tracking"],
           "fps_median": round(q["fps_median"], 1), "fps_mean": round(q["fps_mean"], 1),
           "max_gap_s": round(q["max_gap_s"], 3), "n_frames": y.size, "hold_s": round(q["hold_s"], 2),
           "n_spikes": q["n_spikes"], "sigma_nm": round(sig, 2),
           "flicker": round(float(np.mean(np.abs(r9) > 4 * sig)), 4),
           "repeat": round(float(cnt.max() / y.size), 4),
           "lo_bp_B": round(float(lv.min()), 1), "hi_bp_B": round(float(lv.max()), 1)}
    reason = []
    if not SIG_MIN_NM <= sig <= SIG_MAX_NM:
        reason.append("sigma")
    if rec["flicker"] > FLICKER_MAX:
        reason.append("flicker")
    if rec["repeat"] > REPEAT_MAX:
        reason.append("clipped")
    if lv.min() < GUARD_BP_B[0] or lv.max() > GUARD_BP_B[1]:
        reason.append("LUT level")
    rec["gate_reason"] = "+".join(reason)
    rec["passes"] = not reason
    return rec


def stage_qc(plats, data) -> pd.DataFrame:
    qc = pd.DataFrame(pool_map(qc_one, list(range(len(plats))), plats))
    pr = data["pairs"]
    pr = pr[pr.a_present & pr.b_present]
    ov = {(r.data_id, r.roi_a) for r in pr.itertuples()} | {(r.data_id, r.roi_b) for r in pr.itertuples()}
    qc["overlap_roi"] = [(d, b) in ov for d, b in zip(qc.data_id, qc.bead)]
    qc.to_csv(OUT / "plateaus_qc.csv", index=False)
    pr.to_csv(OUT / "overlapping_roi_pairs.csv", index=False)
    g = qc[qc.passes]
    print(f"qc: {len(qc)} plateaus, {len(g)} pass; user96 {int((g.set == 'user96').sum())} plateaus "
          f"({g[g.set == 'user96'].groupby(['data_id', 'bead']).ngroups} tethers); fail: "
          f"{qc[~qc.passes].gate_reason.value_counts().to_dict()}", flush=True)
    return qc


# ----------------------------------------------------------------- stage settle
SETTLE_BINS = np.arange(0.0, 3.0001, 0.05)


def settle_one(i: int) -> dict:
    """Change points of the OP fit on the WHOLE plateau (t >= 0) vs time since arrival;
    extension relative to the 1.5-3 s median (first 3 s)."""
    q = _P[i]
    y, t, f = q["y"], q["t"], q["force_pn"]
    sig = SB.sigma_diff(y)
    st = fit_levels(t, y, sig, LAM)
    stp = step_table(st, t, f)
    big = stp.cls.isin(["big_fwd", "big_back"]).to_numpy()
    h_all, _ = np.histogram(stp.t_s, SETTLE_BINS)
    h_small, _ = np.histogram(stp.t_s[~big], SETTLE_BINS)
    ref = np.median(y[(t >= 1.5) & (t < 3.0)]) if np.sum((t >= 1.5) & (t < 3.0)) > 5 else np.nan
    idx = np.digitize(t, SETTLE_BINS) - 1
    dev = np.full(SETTLE_BINS.size - 1, np.nan)
    for k in range(SETTLE_BINS.size - 1):
        m = idx == k
        if m.any():
            dev[k] = float(np.mean(y[m] - ref))
    first_small = stp.t_s[~big].min() if (~big).any() else np.nan
    return {"i": i, "force_pn": f, "h_all": h_all, "h_small": h_small, "dev": dev,
            "y_first_bp": float(luz_b(st.level_nm.iloc[0], f)[0])}


def stage_settle(plats, qc) -> pd.DataFrame:
    g = qc[qc.passes]
    res = pool_map(settle_one, g.i.tolist(), plats)
    rows = []
    centers = 0.5 * (SETTLE_BINS[:-1] + SETTLE_BINS[1:])
    for grp, sel in (("10 pN", lambda r: r["force_pn"] == 10.0), (">=18 pN", lambda r: r["force_pn"] >= 18)):
        rr = [r for r in res if sel(r)]
        n = len(rr)
        ha = np.sum([r["h_all"] for r in rr], 0)
        hs = np.sum([r["h_small"] for r in rr], 0)
        dev = np.nanmean(np.vstack([r["dev"] for r in rr]), 0)
        for c, a_, s_, d_ in zip(centers, ha, hs, dev):
            rows.append({"group": grp, "t_mid_s": round(c, 3), "n_plateaus": n,
                         "change_points_per_plateau_per_s": a_ / n / 0.05,
                         "small_change_points_per_plateau_per_s": s_ / n / 0.05,
                         "mean_ext_minus_1p5_3s_median_nm": d_})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "settle_changepoint_density.csv", index=False)
    return out


# ----------------------------------------------------------------- stage detect
EXP_DT = 0.01
DT_EDGES = np.round(np.arange(-3.0, 3.0001, 0.1), 3)
LUZ_EDGES = np.arange(-150.0, 1150.0001, 10.0)


def exposure(t, st, big_t, f):
    """Searchable time on a 0.01 s grid: OP level (Luz_B), open-fork flag, signed
    time to the nearest big event and its kind."""
    tg = np.arange(T_SKIP, t[-1], EXP_DT) + 0.5 * EXP_DT
    edges = np.concatenate([st.t0.to_numpy(), [np.inf]])
    seg = np.clip(np.searchsorted(edges, tg, side="right") - 1, 0, len(st) - 1)
    lv = luz_b(st.level_nm.to_numpy()[seg], f)
    d = tg[:, None] - big_t[None, :]
    j = np.argmin(np.abs(d), axis=1)
    dn = d[np.arange(tg.size), j]
    kind_arr = np.where(big_t[j] == 0.0, 0, 1)            # 0 = arrival, 1 = big jump
    return tg, lv, dn, kind_arr


def window_class(dn, kind):
    """'after arrival' / 'after big' (0-1 s after) / 'before big' (0-1 s before) / '1-2 s' / 'far' (> 2 s)."""
    dn = np.asarray(dn, float)
    kind = np.asarray(kind)
    out = np.full(dn.shape, "1-2 s", dtype=object)
    out[np.abs(dn) > FAR_S] = "far"
    near = np.abs(dn) <= WIN_S
    out[near & (dn >= 0) & (kind == 0)] = "after arrival"
    out[near & (dn >= 0) & (kind == 1)] = "after big"
    out[near & (dn < 0)] = "before big"
    return out


WCLASSES = ("after arrival", "after big", "before big", "1-2 s", "far")


def exposure_summary(tg, lv, dn, kind) -> dict:
    op = (lv >= OPEN_LO) & (lv < OPEN_HI)
    wc = window_class(dn, kind)
    ex = {"exp_total_s": tg.size * EXP_DT, "exp_open_s": op.sum() * EXP_DT,
          "exp_closed_s": (lv < OPEN_LO).sum() * EXP_DT, "exp_full_s": (lv >= OPEN_HI).sum() * EXP_DT}
    for w in WCLASSES:
        ex[f"exp_open_{w}_s"] = float(np.sum(op & (wc == w)) * EXP_DT)
    h_arr, _ = np.histogram(dn[op & (kind == 0)], DT_EDGES)
    h_big, _ = np.histogram(dn[op & (kind == 1)], DT_EDGES)
    occ, _ = np.histogram(lv[op], LUZ_EDGES)
    return ex, h_arr * EXP_DT, h_big * EXP_DT, occ * EXP_DT


def detect_one(i: int):
    q = _P[i]
    f = q["force_pn"]
    t, y = search(q)
    sig = SB.sigma_diff(y)
    st = fit_levels(t, y, sig, LAM)
    stp = step_table(st, t, f)
    big_t = big_events(stp)
    stp = relate(stp, big_t, st, f)
    db, kpl, nw = [], [], []
    for r in stp.itertuples():
        a, b, c = event_dbic(t, y, st, int(r.k), sig) if r.cls.startswith("small") else (np.nan, np.nan, 0)
        db.append(a)
        kpl.append(b)
        nw.append(c)
    stp["dbic_pl_minus_step"] = db
    stp["pl_knots"] = kpl
    stp["n_window"] = nw
    stp["t_abs_s"] = q["t_file_s"] + t[stp.idx.to_numpy(int)] if len(stp) else []
    stp["tol_s"] = [max(t[min(k + 1, t.size - 1)] - t[k], t[k] - t[max(k - 1, 0)]) + 1e-6
                    for k in stp.idx.to_numpy(int)]
    tg, lv, dn, kind = exposure(t, st, big_t, f)
    ex, h_arr, h_big, occ = exposure_summary(tg, lv, dn, kind)
    lv_st = luz_b(st.level_nm.to_numpy(), f)
    base = {"i": i, "data_id": q["data_id"], "bead": q["bead"], "step": q["step"], "force_pn": f,
            "set": "user96" if q["user96"] else "other_clean", "sigma_nm": round(sig, 2),
            "fps_mean": round(q["fps_mean"], 1), "fps_median": round(q["fps_median"], 1)}
    if len(stp):
        stp["wclass"] = window_class(stp.dt_nearest_big_s.to_numpy(),
                                     (stp.nearest_big_kind == "big jump").astype(int).to_numpy())
        stp["open_fork"] = (stp.from_bp >= OPEN_LO) & (stp.from_bp < OPEN_HI)
    pl = {**base, "n_levels": len(st), "first_level_bp": round(float(lv_st[0]), 1),
          "max_level_bp": round(float(lv_st.max()), 1), "last_level_bp": round(float(lv_st[-1]), 1),
          "n_big_fwd": int((stp.cls == "big_fwd").sum()) if len(stp) else 0,
          "n_big_back": int((stp.cls == "big_back").sum()) if len(stp) else 0,
          "opened": bool(lv_st.max() >= OPEN_LO), **ex}
    return [{**base, **r} for r in stp.to_dict("records")], pl, h_arr, h_big, occ


def stage_detect(plats, qc):
    g = qc[qc.passes]
    res = pool_map(detect_one, g.i.tolist(), plats, chunksize=2)
    steps = pd.DataFrame([r for s_, _, _, _, _ in res for r in s_])
    pls = pd.DataFrame([p_ for _, p_, _, _, _ in res])
    expo = {"i": np.array([p_["i"] for _, p_, _, _, _ in res]),
            "dt_arrival": np.vstack([a for _, _, a, _, _ in res]),
            "dt_big": np.vstack([b for _, _, _, b, _ in res]),
            "occupancy": np.vstack([o for _, _, _, _, o in res])}
    steps.to_csv(OUT / "steps_all.csv.gz", index=False)
    pls.to_csv(OUT / "plateau_summary.csv", index=False)
    np.savez_compressed(OUT / "_exposure.npz", **expo, dt_edges=DT_EDGES, luz_edges=LUZ_EDGES)
    print(f"detect: {len(steps)} change points on {len(pls)} plateaus; classes "
          f"{steps.cls.value_counts().to_dict()}", flush=True)
    return steps, pls, expo


# ----------------------------------------------------------------- stage nulls (N2 / N3-PL)
def skeletons(t, y, f, sig, st, stp):
    """N2: flat level per stretch between |step| >= 100 bp B; N3: continuous PL
    (0-3 knots on a >= KNOT_GRID_S grid) per stretch; residuals of the lambda = 3 fit."""
    n = y.size
    cut = [0] + sorted(stp.loc[stp.cls.isin(["big_fwd", "big_back"]), "idx"].astype(int).tolist()) + [n]
    n2 = np.empty(n)
    n3 = np.empty(n)
    knots_used = []
    for a, b in zip(cut[:-1], cut[1:]):
        ys, ts = y[a:b], t[a:b]
        n2[a:b] = ys.mean()
        if b - a < 8:
            n3[a:b] = ys.mean()
            continue
        span = ts[-1] - ts[0]
        step_ = max(KNOT_GRID_S, span / 60.0)
        cand = np.arange(ts[0] + step_, ts[-1] - 0.5 * step_, step_)
        if cand.size == 0:
            n3[a:b] = np.polyval(np.polyfit(ts, ys, 1), ts)
            continue
        _, k, _, kn, curve = pl_fit(ts, ys, cand, sig, want_curve=True)
        n3[a:b] = curve
        knots_used.append(k)
    ref = fit_levels(t, y, sig, LAM_REF)
    res = y - fit_curve(ref, n)
    return n2, n3, res, knots_used


def nulls_one(i: int):
    q = _P[i]
    f = q["force_pn"]
    t, y = search(q)
    sig = SB.sigma_diff(y)
    st = fit_levels(t, y, sig, LAM)
    stp = step_table(st, t, f)
    big_t = big_events(stp)
    n2, n3, res, kn = skeletons(t, y, f, sig, st, stp)
    rng = np.random.default_rng(91_000 + i)
    rows, cnt = [], []
    for name, base in (("N2", n2), ("N3", n3)):
        for rep in range(N_NULL):
            ys = base + SB.block_boot(res, y.size, rng)
            s_ = SB.sigma_diff(ys)
            sst = fit_levels(t, ys, s_, LAM)
            sp = step_table(sst, t, f)
            small = sp[sp.cls.str.startswith("small")]
            for r in small.itertuples():
                d = r.t_s - big_t
                j = int(np.argmin(np.abs(d)))
                db, kpl, nw = event_dbic(t, ys, sst, int(r.k), s_)
                rows.append({"i": i, "force_pn": f, "set": "user96" if q["user96"] else "other_clean",
                             "null": name, "rep": rep, "t_s": r.t_s, "step_bp": r.step_bp, "se_bp": r.se_bp,
                             "from_bp": r.from_bp, "to_bp": r.to_bp, "cls": r.cls,
                             "dwell_before_s": r.dwell_before_s, "dwell_after_s": r.dwell_after_s,
                             "dt_nearest_big_s": float(d[j]), "nearest_big_kind": "arrival" if big_t[j] == 0 else "big jump",
                             "dbic_pl_minus_step": db, "pl_knots": kpl})
            op = (small.from_bp >= OPEN_LO) & (small.from_bp < OPEN_HI)
            cnt.append({"i": i, "force_pn": f, "null": name, "rep": rep,
                        "n_small": len(small), "n_small_open": int(op.sum()),
                        "n_small_fwd_open": int((op & (small.cls == "small_fwd")).sum()),
                        "n_small_back_open": int((op & (small.cls == "small_back")).sum()),
                        "n_big_fwd": int((sp.cls == "big_fwd").sum())})
    return rows, cnt, kn


def stage_nulls(plats, qc):
    g = qc[qc.passes & (qc.force_pn >= 18)]
    res = pool_map(nulls_one, g.i.tolist(), plats, chunksize=2)
    ev = pd.DataFrame([r for rr, _, _ in res for r in rr])
    cnt = pd.DataFrame([c for _, cc, _ in res for c in cc])
    kn = np.concatenate([np.asarray(k, int) for _, _, k in res if len(k)])
    ev.to_csv(OUT / "null_small_jumps.csv.gz", index=False)
    cnt.to_csv(OUT / "null_counts.csv.gz", index=False)
    print(f"nulls: {len(ev)} small jumps on {len(g)} plateaus x {N_NULL} reps; N3 knots used "
          f"{np.bincount(kn, minlength=4).tolist()}", flush=True)
    return ev, cnt


# ----------------------------------------------------------------- stage inject
def inject_one(i: int):
    """Abrupt single steps of A bp B (forward) inserted at a random time into this
    plateau's own N2 base (own clock, own bootstrap noise); detected if a change
    point lies within 2 frames with a size within +-50 %."""
    q = _P[i]
    f = q["force_pn"]
    g = gain_b(f)
    t, y = search(q)
    sig = SB.sigma_diff(y)
    st = fit_levels(t, y, sig, LAM)
    stp = step_table(st, t, f)
    big_t = big_events(stp)
    n = y.size                                   # N2 base only (no PL fit needed here)
    cut = [0] + sorted(stp.loc[stp.cls.isin(["big_fwd", "big_back"]), "idx"].astype(int).tolist()) + [n]
    base = np.empty(n)
    for a, b in zip(cut[:-1], cut[1:]):
        base[a:b] = y[a:b].mean()
    ref = fit_levels(t, y, sig, LAM_REF)
    res = y - fit_curve(ref, n)
    rng = np.random.default_rng(55_000 + i)
    cuts_t = np.array([t[c] for c in cut[1:-1]])
    rows = []
    for amp in INJ_AMPS_BP:
        for rep in range(INJ_REPS):
            ts_ = rng.uniform(T_SKIP + 0.05, t[-1] - 0.05)
            s_idx = int(np.searchsorted(t, ts_))
            if s_idx < 1 or s_idx >= n:
                continue
            ys = base + SB.block_boot(res, n, rng)
            ys[s_idx:] += amp * g
            s_ = SB.sigma_diff(ys)
            sst = fit_levels(t, ys, s_, LAM)
            sp = step_table(sst, t, f)
            hit = sp[(np.abs(sp.idx - s_idx) <= 2) & (np.abs(sp.step_nm - amp * g) <= 0.5 * amp * g)]
            prev_c = np.concatenate([[T_SKIP], cuts_t[cuts_t < ts_]]).max()
            next_c = np.concatenate([[t[-1]], cuts_t[cuts_t > ts_]]).min()
            d = ts_ - big_t
            j = int(np.argmin(np.abs(d)))
            row = {"i": i, "force_pn": f, "set": "user96" if q["user96"] else "other_clean",
                   "sigma_nm": round(sig, 2), "fps_mean": round(q["fps_mean"], 1), "amp_bp": amp,
                   "amp_nm": amp * g, "t_inj_s": ts_, "dwell_before_s": ts_ - prev_c, "dwell_after_s": next_c - ts_,
                   "dt_nearest_big_s": float(d[j]), "nearest_big_kind": "arrival" if big_t[j] == 0 else "big jump",
                   "detected": len(hit) > 0}
            if len(hit):
                r = hit.iloc[0]
                db, kpl, nw = event_dbic(t, ys, sst, int(r.k), s_)
                row.update(est_bp=r.step_bp, est_cls=r.cls, dbic_pl_minus_step=db, pl_knots=kpl)
            rows.append(row)
    return rows


def stage_inject(plats, qc):
    g = qc[qc.passes & (qc.force_pn >= 18)]
    res = pool_map(inject_one, g.i.tolist(), plats, chunksize=2)
    inj = pd.DataFrame([r for rr in res for r in rr])
    inj.to_csv(OUT / "injection_recall.csv.gz", index=False)
    print(f"inject: {len(inj)} injections, recall {inj.groupby('amp_bp').detected.mean().round(3).to_dict()}",
          flush=True)
    return inj




# ----------------------------------------------------------------- stage injreal (real-noise calibration)
REAL_AMPS_BP = (20.0, 40.0, 60.0, 80.0)
REAL_AMPS_NM_10PN = (20.0, 30.0, 40.0, 60.0, 80.0)
SHAPES = (("abrupt", 0.0), ("ramp 0.3 s", 0.3), ("ramp 1.0 s", 1.0))
REAL_REPS = 2


def injreal_one(i: int):
    """Calibration in REAL noise: an abrupt step, or a linear ramp of 0.3 / 1.0 s,
    of A bp B added to the real trace (its own drift, low-frequency noise, own
    clock and its own fork events) at a random time. Every change point found
    inside the event (+-2 frames) with 25-150 % of A (abrupt: 50-150 %) is a
    detected piece; its dBIC is computed exactly as for the real small jumps.
    10 pN plateaus (closed fork): abrupt steps of 20-80 nm, recall only."""
    q = _P[i]
    f = q["force_pn"]
    g = gain_b(f)
    t, y = search(q)
    sig0 = SB.sigma_diff(y)
    st0 = fit_levels(t, y, sig0, LAM)
    stp0 = step_table(st0, t, f)
    big_t = big_events(stp0)
    rng = np.random.default_rng(77_000 + i)
    n = y.size
    rows = []
    if f == 10.0:
        jobs = [(a, "abrupt", 0.0) for a in REAL_AMPS_NM_10PN for _ in range(3)]
    else:
        jobs = [(a * g, sh, tau) for a in REAL_AMPS_BP for sh, tau in SHAPES for _ in range(REAL_REPS)]
    for amp_nm, shape, tau in jobs:
        ts_ = rng.uniform(T_SKIP + 0.3, max(T_SKIP + 0.31, t[-1] - tau - 0.3))
        add = amp_nm * (np.clip((t - ts_) / tau, 0.0, 1.0) if tau > 0 else (t > ts_).astype(float))
        ys = y + add
        s_ = SB.sigma_diff(ys)
        sst = fit_levels(t, ys, s_, LAM)
        sp = step_table(sst, t, f)
        lo_i = int(np.searchsorted(t, ts_)) - 2
        hi_i = int(np.searchsorted(t, ts_ + tau)) + 2
        frac_lo = 0.5 if tau == 0 else 0.25
        hit = sp[(sp.idx >= lo_i) & (sp.idx <= hi_i) & (sp.step_nm >= frac_lo * amp_nm) & (sp.step_nm <= 1.5 * amp_nm)]
        d = ts_ - big_t
        j = int(np.argmin(np.abs(d)))
        base = {"i": i, "force_pn": f, "set": "user96" if q["user96"] else "other_clean", "shape": shape,
                "amp_nm": amp_nm, "amp_bp": amp_nm / g, "sigma_nm": round(sig0, 2), "t_inj_s": ts_,
                "dt_nearest_big_s": float(d[j]), "detected": len(hit) > 0, "n_pieces": len(hit)}
        if not len(hit):
            rows.append(base)
            continue
        for r in hit.itertuples():
            db, kpl, nw = event_dbic(t, ys, sst, int(r.k), s_) if abs(r.step_bp) < BIG_BP else (np.nan, np.nan, 0)
            rows.append({**base, "piece_bp": r.step_bp, "piece_nm": r.step_nm, "piece_cls": r.cls,
                         "dbic_pl_minus_step": db, "pl_knots": kpl})
    return rows


def stage_injreal(plats, qc):
    g = qc[qc.passes]
    res = pool_map(injreal_one, g.i.tolist(), plats, chunksize=2)
    inj = pd.DataFrame([r for rr in res for r in rr])
    inj.to_csv(OUT / "injection_real_noise.csv.gz", index=False)
    ev = inj.drop_duplicates(["i", "shape", "amp_nm", "t_inj_s"])
    print("injreal recall:", ev.groupby(["shape", ev.force_pn == 10.0]).detected.mean().round(3).to_dict(), flush=True)
    return inj


# ----------------------------------------------------------------- stage stats
FORCES = (18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9, 46.0)
N_BOOT = 2000


def tether_key(df):
    return df.data_id.astype(int).astype(str) + "_" + df.bead.astype(str)


def boot_ratio(num, den, teth, rng, n=N_BOOT):
    """Ratio sum(num) / sum(den) with a 95 % CI from resampling tethers."""
    num, den = np.asarray(num, float), np.asarray(den, float)
    u, inv = np.unique(np.asarray(teth), return_inverse=True)
    if u.size == 0 or den.sum() <= 0:
        return np.nan, np.nan, np.nan
    nu = np.bincount(inv, num, u.size)
    de = np.bincount(inv, den, u.size)
    pick = rng.integers(0, u.size, (n, u.size))
    with np.errstate(invalid="ignore", divide="ignore"):
        b = nu[pick].sum(1) / de[pick].sum(1)
    return float(nu.sum() / de.sum()), float(np.nanpercentile(b, 2.5)), float(np.nanpercentile(b, 97.5))


def ungated_changes(plats, idx_list):
    """OP (lambda = 4) change points of plateaus, gated or not (for the overlap check)."""
    out = {}
    for i in idx_list:
        q = plats[i]
        t, y = search(q)
        if y.size < 2 * MMIN + 2:
            continue
        sig = SB.sigma_diff(y)
        st = fit_levels(t, y, sig, LAM)
        out[i] = q["t_file_s"] + t[st.a.to_numpy()[1:].astype(int)]
    return out


def coincidence(steps, plats, data, qc, rng):
    """Copy check. (1) overlapping ROI pairs (< 60 px, both beads present): a small
    jump of A at the same frame (+-1) as any change point of B (B gated or not) is a
    copy suspect. (2) all other beads of the same file / step (gated): same-frame
    coincidence vs a time-shift null (global artefacts)."""
    pairs = data["pairs"]
    pairs = pairs[pairs.a_present & pairs.b_present]
    key = {(q["data_id"], q["bead"], q["step"]): i for i, q in enumerate(plats)}
    nb = {}
    for r in pairs.itertuples():
        nb.setdefault((r.data_id, r.roi_a), set()).add(r.roi_b)
        nb.setdefault((r.data_id, r.roi_b), set()).add(r.roi_a)
    need = set()
    sm = steps[steps.cls.str.startswith("small") & (steps.force_pn >= 18)]
    for r in sm.itertuples():
        for b in nb.get((r.data_id, r.bead), ()):
            j = key.get((r.data_id, b, r.step))
            if j is not None:
                need.add(j)
    ch = ungated_changes(plats, sorted(need))
    copy_flag, copy_shift = [], []
    for r in sm.itertuples():
        hit = False
        hit_s = False
        sh = rng.choice([-1, 1]) * rng.uniform(0.5, 3.0)
        for b in nb.get((r.data_id, r.bead), ()):
            j = key.get((r.data_id, b, r.step))
            if j is None or j not in ch:
                continue
            hit |= bool(np.any(np.abs(ch[j] - r.t_abs_s) <= r.tol_s))
            hit_s |= bool(np.any(np.abs(ch[j] - (r.t_abs_s + sh)) <= r.tol_s))
        copy_flag.append(hit)
        copy_shift.append(hit_s)
    steps = steps.copy()
    steps["copy_suspect"] = False
    steps.loc[sm.index, "copy_suspect"] = copy_flag
    steps["has_overlap_neighbour"] = [(d, b) in nb for d, b in zip(steps.data_id, steps.bead)]
    # (2) global: any other gated bead of the file, same step, same frame
    allc = steps[steps.force_pn >= 18][["data_id", "bead", "step", "t_abs_s"]]
    glob, glob_null = [], []
    grp = {k: g for k, g in allc.groupby(["data_id", "step"])}
    for r in sm.itertuples():
        g = grp[(r.data_id, r.step)]
        o = g[(g.bead != r.bead) & ~g.bead.isin(nb.get((r.data_id, r.bead), set()))].t_abs_s.to_numpy()
        glob.append(int(np.sum(np.abs(o - r.t_abs_s) <= r.tol_s)))
        sh = rng.choice([-1, 1]) * rng.uniform(0.5, 3.0)
        glob_null.append(int(np.sum(np.abs(o - (r.t_abs_s + sh)) <= r.tol_s)))
    steps["n_other_beads_same_frame"] = np.nan
    steps.loc[sm.index, "n_other_beads_same_frame"] = glob
    summ = {"small_jumps_checked": int(len(sm)),
            "in_overlap_pairs": int(sm.apply(lambda r: (r.data_id, r.bead) in nb, axis=1).sum()),
            "copy_suspects": int(np.sum(copy_flag)), "copy_suspects_shift_null": int(np.sum(copy_shift)),
            "coincident_with_other_bead": int(np.sum(np.array(glob) > 0)),
            "coincident_with_other_bead_shift_null": int(np.sum(np.array(glob_null) > 0)),
            "coincident_with_2plus_other_beads": int(np.sum(np.array(glob) > 1)),
            "coincident_with_2plus_shift_null": int(np.sum(np.array(glob_null) > 1))}
    return steps, summ


def recall_table(inj: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (s_, f), g in inj.groupby(["set", "force_pn"]):
        for a, h in g.groupby("amp_bp"):
            longd = h[(h.dwell_before_s >= 0.5) & (h.dwell_after_s >= 0.5)]
            rows.append({"set": s_, "force_pn": f, "amp_bp": a, "amp_nm": round(float(h.amp_nm.mean()), 1),
                         "n_inj": len(h), "recall": round(float(h.detected.mean()), 3),
                         "recall_both_dwells_ge_0p5s": round(float(longd.detected.mean()), 3) if len(longd) else np.nan,
                         "recall_min_dwell_lt_0p3s": round(float(h[np.minimum(h.dwell_before_s, h.dwell_after_s) < 0.3]
                                                               .detected.mean()), 3),
                         "recall_within_1s_of_big": round(float(h[np.abs(h.dt_nearest_big_s) <= 1].detected.mean()), 3),
                         "recall_far_gt_2s": round(float(h[np.abs(h.dt_nearest_big_s) > 2].detected.mean()), 3),
                         "median_sigma_nm": round(float(h.sigma_nm.median()), 1),
                         "median_fps_mean": round(float(h.fps_mean.median()), 1)})
    return pd.DataFrame(rows)


def a80(rt: pd.DataFrame, col="recall", level=0.8) -> float:
    """Smallest amplitude (bp B, interpolated on the injection grid) from which recall stays >= level."""
    r = rt.sort_values("amp_bp")
    a, v = r.amp_bp.to_numpy(), r[col].to_numpy()
    for j in range(a.size):
        if np.all(v[j:] >= level):
            if j == 0:
                return float(a[0])
            return float(np.interp(level, [v[j - 1], v[j]], [a[j - 1], a[j]])) if v[j] > v[j - 1] else float(a[j])
    return np.nan


def recall_at(rt: pd.DataFrame, amp):
    r = rt.sort_values("amp_bp")
    return np.interp(np.abs(amp), r.amp_bp.to_numpy(), r.recall.to_numpy())


def group_force(df, f):
    if f == ">=18":
        return df[df.force_pn >= 18]
    return df[np.isclose(df.force_pn, float(f))]


def stage_stats(plats, qc, data):
    rng = np.random.default_rng(2026)
    steps = pd.read_csv(OUT / "steps_all.csv.gz")
    pls = pd.read_csv(OUT / "plateau_summary.csv")
    nev = pd.read_csv(OUT / "null_small_jumps.csv.gz")
    ncnt = pd.read_csv(OUT / "null_counts.csv.gz")
    inj = pd.read_csv(OUT / "injection_recall.csv.gz")
    injr = pd.read_csv(OUT / "injection_real_noise.csv.gz")
    expo = np.load(OUT / "_exposure.npz")
    stats = {"params": {"lambda": LAM, "T_skip_s": T_SKIP, "big_bp_B": BIG_BP, "open_fork_bp_B": [OPEN_LO, OPEN_HI],
                        "window_s": WIN_S, "far_s": FAR_S, "n_null": N_NULL, "knot_grid_s": KNOT_GRID_S,
                        "disc_window_s": DISC_WIN_S, "guard_bp_B": GUARD_BP_B}}
    # ---- gate summary
    g = qc[qc.passes]
    stats["gate"] = {s_: {"plateaus_collected": int((qc.set == s_).sum()), "plateaus_pass": int((g.set == s_).sum()),
                          "plateaus_pass_ge18": int(((g.set == s_) & (g.force_pn >= 18)).sum()),
                          "tethers_pass": int(g[g.set == s_].groupby(["data_id", "bead"]).ngroups),
                          "fail_reasons": qc[(qc.set == s_) & ~qc.passes].gate_reason.value_counts().to_dict(),
                          "median_sigma_nm": round(float(g[g.set == s_].sigma_nm.median()), 1)} for s_ in SETS}
    # ---- copy / coincidence checks
    steps, coin = coincidence(steps, plats, data, qc, rng)
    stats["coincidence"] = coin
    steps["tether"] = tether_key(steps)
    pls["tether"] = tether_key(pls)
    steps["use"] = ~steps.copy_suspect.astype(bool)
    steps.to_csv(OUT / "steps_all_flagged.csv.gz", index=False)
    pl18 = pls[pls.force_pn >= 18].copy()
    global OVERLAP_SET
    pr_ = data["pairs"]
    pr_ = pr_[pr_.a_present & pr_.b_present]
    OVERLAP_SET = {(r.data_id, r.roi_a) for r in pr_.itertuples()} | {(r.data_id, r.roi_b) for r in pr_.itertuples()}
    sm = steps[steps.cls.str.startswith("small") & (steps.force_pn >= 18) & steps.use]
    smo = sm[sm.open_fork.astype(bool)]
    # ---- detection limits
    rt = recall_table(inj)
    rt.to_csv(OUT / "detection_recall_by_force.csv", index=False)
    lim = []
    for (s_, f), r in rt.groupby(["set", "force_pn"]):
        lim.append({"set": s_, "force_pn": f, "A80_bp_B": a80(r), "A50_bp_B": a80(r, level=0.5),
                    "A80_bp_B_long_dwells": a80(r, "recall_both_dwells_ge_0p5s"),
                    "A80_nm": a80(r) * gain_b(f) if np.isfinite(a80(r)) else np.nan,
                    "median_sigma_nm": float(r.median_sigma_nm.iloc[0]), "median_fps_mean": float(r.median_fps_mean.iloc[0])})
    lim = pd.DataFrame(lim)
    lim.to_csv(OUT / "detection_limits_by_force.csv", index=False)
    # ---- per-force probabilities and rates
    ncnt = ncnt.merge(pls[["i", "set", "tether", "opened"]], on="i", how="left")
    rows = []
    for s_ in SETS:
        for f in FORCES + (">=18",):
            P_ = group_force(pl18[pl18.set == s_], f)
            P_ = P_[P_.opened]
            if not len(P_):
                continue
            E_ = smo[smo.i.isin(P_.i)]
            n_f = E_[E_.cls == "small_fwd"].groupby("i").size().reindex(P_.i, fill_value=0).to_numpy()
            n_b = E_[E_.cls == "small_back"].groupby("i").size().reindex(P_.i, fill_value=0).to_numpy()
            n_all = n_f + n_b
            nc = ncnt[ncnt.i.isin(P_.i)]
            null_p = {}
            null_rate = {}
            for nn in ("N2", "N3"):
                c = nc[nc.null == nn]
                null_p[nn] = float((c.n_small_open > 0).mean()) if len(c) else np.nan
                null_rate[nn] = {k: float(c[col].sum() / N_NULL / P_.exp_open_s.sum())
                                 for k, col in (("all", "n_small_open"), ("fwd", "n_small_fwd_open"),
                                                ("back", "n_small_back_open"))} if len(c) else {}
            teth = P_.tether.to_numpy()
            p_obs, p_lo, p_hi = boot_ratio(n_all > 0, np.ones(len(P_)), teth, rng)
            r_all = boot_ratio(n_all, P_.exp_open_s, teth, rng)
            r_f = boot_ratio(n_f, P_.exp_open_s, teth, rng)
            r_b = boot_ratio(n_b, P_.exp_open_s, teth, rng)
            pn = null_p["N2"]
            # detectable range: |step| >= A80 at this force (recall-corrected, N2-subtracted)
            fl = 30.0 if f == ">=18" else float(f)
            rts = rt[(rt.set == s_) & np.isclose(rt.force_pn, fl)] if f != ">=18" else None
            row = {"set": s_, "force_pn": f, "n_tethers": int(np.unique(teth).size), "n_opened_plateaus": len(P_),
                   "open_fork_s": round(float(P_.exp_open_s.sum()), 1),
                   "n_small_fwd": int(n_f.sum()), "n_small_back": int(n_b.sum()),
                   "P_plateau_ge1_small": round(p_obs, 3), "P_ci_lo": round(p_lo, 3), "P_ci_hi": round(p_hi, 3),
                   "P_N2_null": round(pn, 3), "P_N3_null": round(null_p["N3"], 3),
                   "P_corrected_N2": round((p_obs - pn) / (1 - pn), 3) if np.isfinite(pn) and pn < 1 else np.nan,
                   "rate_small_per_s": round(r_all[0], 4), "rate_ci_lo": round(r_all[1], 4), "rate_ci_hi": round(r_all[2], 4),
                   "rate_fwd_per_s": round(r_f[0], 4), "rate_back_per_s": round(r_b[0], 4),
                   "rate_N2_per_s": round(null_rate["N2"].get("all", np.nan), 4),
                   "rate_N2_fwd_per_s": round(null_rate["N2"].get("fwd", np.nan), 4),
                   "rate_N2_back_per_s": round(null_rate["N2"].get("back", np.nan), 4),
                   "rate_N3_per_s": round(null_rate["N3"].get("all", np.nan), 4),
                   "rate_net_N2_per_s": round(r_all[0] - null_rate["N2"].get("all", np.nan), 4),
                   "fwd_back_ratio": round(n_f.sum() / max(n_b.sum(), 1), 2)}
            if rts is not None and len(rts):
                A = a80(rts)
                row["A80_bp_B"] = round(A, 1) if np.isfinite(A) else np.nan
                if np.isfinite(A):
                    inr = E_[E_.step_bp.abs() >= A]
                    w = 1.0 / np.maximum(recall_at(rts, inr.step_bp.abs()), 0.05)
                    nn2 = nev[(nev.null == "N2") & nev.i.isin(P_.i) & (nev.step_bp.abs() >= A)
                              & (nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI)]
                    w2 = 1.0 / np.maximum(recall_at(rts, nn2.step_bp.abs()), 0.05)
                    nd = inr.groupby("i").size().reindex(P_.i, fill_value=0).to_numpy()
                    pd_, pdl, pdh = boot_ratio(nd > 0, np.ones(len(P_)), teth, rng)
                    row.update(P_plateau_ge1_detectable=round(pd_, 3), P_det_ci=f"{pdl:.3f}-{pdh:.3f}")
                    row.update(n_small_detectable=len(inr),
                               rate_detectable_raw_per_s=round(len(inr) / P_.exp_open_s.sum(), 4),
                               rate_detectable_recall_corr_per_s=round(float(w.sum()) / P_.exp_open_s.sum(), 4),
                               rate_detectable_N2_recall_corr_per_s=round(float(w2.sum()) / N_NULL / P_.exp_open_s.sum(), 4))
                    row["rate_detectable_net_per_s"] = round(row["rate_detectable_recall_corr_per_s"]
                                                             - row["rate_detectable_N2_recall_corr_per_s"], 4)
            rows.append(row)
    pf = pd.DataFrame(rows)
    for s_ in SETS:                       # pooled >= 18 pN detectable-range rates: sum over forces
        q = pf[(pf.set == s_) & (pf.force_pn.astype(str) != ">=18") & pf.A80_bp_B.notna()]
        T = q.open_fork_s.sum()
        k = pf.index[(pf.set == s_) & (pf.force_pn.astype(str) == ">=18")]
        for c in ("rate_detectable_raw_per_s", "rate_detectable_recall_corr_per_s", "rate_detectable_N2_recall_corr_per_s",
                  "rate_detectable_net_per_s"):
            pf.loc[k, c] = round(float((q[c] * q.open_fork_s).sum() / T), 4)
        pf.loc[k, "n_small_detectable"] = q.n_small_detectable.sum()
        pf.loc[k, "P_plateau_ge1_detectable"] = round(float(np.average(q.P_plateau_ge1_detectable, weights=q.n_opened_plateaus)), 3)
        pf.loc[k, "A80_bp_B"] = round(float(np.average(q.A80_bp_B, weights=q.open_fork_s)), 1)
    pf.to_csv(OUT / "small_jump_probability_by_force.csv", index=False)
    # ---- controls: closed fork, fully open, 10 pN
    ctrl = []
    for s_ in SETS:
        P_ = pls[(pls.set == s_)]
        E_ = steps[(steps.set == s_) & steps.cls.str.startswith("small") & steps.use]
        for lab, pp, ee, ex in (
                ("10 pN plateaus (closed fork)", P_[P_.force_pn == 10.0], E_[E_.force_pn == 10.0], "exp_total_s"),
                (">=18 pN, level < 50 bp B (closed)", P_[P_.force_pn >= 18], E_[(E_.force_pn >= 18) & (E_.from_bp < OPEN_LO)], "exp_closed_s"),
                (">=18 pN, level >= 900 bp B (full)", P_[P_.force_pn >= 18], E_[(E_.force_pn >= 18) & (E_.from_bp >= OPEN_HI)], "exp_full_s"),
                (">=18 pN, open fork 50-900 bp B", P_[P_.force_pn >= 18], E_[(E_.force_pn >= 18) & E_.open_fork.astype(bool)], "exp_open_s")):
            T = float(pp[ex].sum())
            ctrl.append({"set": s_, "control": lab, "time_s": round(T, 1), "n_small": len(ee),
                         "rate_per_s": round(len(ee) / T, 4) if T > 0 else np.nan,
                         "median_sigma_nm": round(float(pp.sigma_nm.median()), 1)})
    ctrl = pd.DataFrame(ctrl)
    ctrl.to_csv(OUT / "controls_no_fork.csv", index=False)
    # ---- per big event: P(small within 1 s before / after), enrichment
    pb = per_big(steps, pls, nev, rng)
    pb.to_csv(OUT / "per_big_event.csv", index=False)
    en = enrichment(smo, pl18, nev, inj, rng)
    en.to_csv(OUT / "enrichment_by_window.csv", index=False)
    # ---- delta-t to nearest big event (rate per s of open-fork time)
    dtt = dt_table(smo, pl18, nev, expo)
    dtt.to_csv(OUT / "rate_vs_dt_to_big.csv", index=False)
    # ---- empirical no-fork controls, amplitude-matched and recall-corrected (nm)
    cev = control_events(steps)
    cmp_ = control_comparison(steps, pls, injr, rng)
    cmp_.to_csv(OUT / "controls_amplitude_matched.csv", index=False)
    # ---- discreteness
    disc, dsum, cal, theta = discreteness(smo, nev, inj, injr, cev, rng)
    dsum.to_csv(OUT / "discreteness_summary.csv", index=False)
    stats["discreteness_theta_by_abs_step_bp"] = {f"{AMP_BINS_DISC[k]:g}-{AMP_BINS_DISC[k + 1]:g}": round(v, 2)
                                                  for k, v in theta.items()}
    stats["discreteness_calibration"] = cal.to_dict(orient="records")
    stats["controls_amplitude_matched"] = cmp_.to_dict(orient="records")
    # ---- sensitivity: drop small jumps that coincide with a change in another bead of the file
    sens = []
    for s_ in SETS:
        P_ = pl18[(pl18.set == s_) & pl18.opened]
        for lab, E_ in (("main", smo[smo.i.isin(P_.i)]),
                        ("no same-frame change in another bead", smo[smo.i.isin(P_.i) & (smo.n_other_beads_same_frame == 0)]),
                        ("tethers without an overlapping neighbour", smo[smo.i.isin(P_.i) & ~smo.has_overlap_neighbour.astype(bool)])):
            PP = P_ if not lab.startswith("tethers") else P_[~P_.tether.isin(smo[smo.has_overlap_neighbour.astype(bool)].tether.unique())
                                                              & ~P_.apply(lambda r: (r.data_id, r.bead) in OVERLAP_SET, axis=1)]
            EE = E_[E_.i.isin(PP.i)]
            n = EE.groupby("i").size().reindex(PP.i, fill_value=0).to_numpy()
            sens.append({"set": s_, "variant": lab, "n_plateaus": len(PP), "n_small": int(n.sum()),
                         "P_plateau_ge1": round(float(np.mean(n > 0)), 3),
                         "rate_per_s": round(float(n.sum() / PP.exp_open_s.sum()), 4)})
    sens = pd.DataFrame(sens)
    sens.to_csv(OUT / "sensitivity_coincidence_overlap.csv", index=False)
    stats["sensitivity"] = sens.to_dict(orient="records")
    # ---- settling window (from stage settle: OP on the whole plateau, both sets)
    se = pd.read_csv(OUT / "settle_changepoint_density.csv")
    hi18 = se[se.group == ">=18 pN"]
    stats["settle"] = {
        "small_cp_per_plateau_per_s_0_0.3s": round(float(hi18[hi18.t_mid_s < 0.3].small_change_points_per_plateau_per_s.mean()), 3),
        "small_cp_per_plateau_per_s_0.3_1s": round(float(hi18[(hi18.t_mid_s > 0.3) & (hi18.t_mid_s < 1)].small_change_points_per_plateau_per_s.mean()), 3),
        "small_cp_per_plateau_per_s_2_3s": round(float(hi18[hi18.t_mid_s > 2].small_change_points_per_plateau_per_s.mean()), 3),
        "ext_dev_10pN_first_0.3s_nm": round(float(se[(se.group == "10 pN") & (se.t_mid_s < 0.3)].mean_ext_minus_1p5_3s_median_nm.abs().max()), 1),
        "ext_dev_ge18_t0_nm": round(float(hi18.mean_ext_minus_1p5_3s_median_nm.iloc[0]), 1),
        "ext_dev_ge18_0.3s_nm": round(float(hi18[np.isclose(hi18.t_mid_s, 0.325)].mean_ext_minus_1p5_3s_median_nm.iloc[0]), 1),
        "ext_dev_ge18_1s_nm": round(float(hi18[np.isclose(hi18.t_mid_s, 1.025)].mean_ext_minus_1p5_3s_median_nm.iloc[0]), 1)}
    # ---- amplitudes / periodicity
    amp, scan = amplitude_stats(smo, nev, rng)
    amp.to_csv(OUT / "amplitude_tests.csv", index=False)
    scan.to_csv(OUT / "periodicity_scan.csv", index=False)
    # ---- positions
    pos, occ_null = position_stats(smo, pl18, expo, inj, rng)
    pos.to_csv(OUT / "position_tests.csv", index=False)
    occ_null.to_csv(OUT / "position_occupancy_null.csv", index=False)
    rel = relative_landing(smo)
    rel.to_csv(OUT / "relative_to_landing_level.csv", index=False)
    smo.to_csv(OUT / "small_jumps_open_fork.csv", index=False)
    # ---- summary
    stats["detection_limits"] = lim.round(1).to_dict(orient="records")
    stats["per_force"] = pf.to_dict(orient="records")
    stats["controls"] = ctrl.to_dict(orient="records")
    stats["per_big_event"] = pb.to_dict(orient="records")
    stats["enrichment"] = en.to_dict(orient="records")
    stats["discreteness"] = dsum.to_dict(orient="records")
    stats["amplitude"] = amp.to_dict(orient="records")
    stats["positions"] = pos.to_dict(orient="records")
    (OUT / "stats.json").write_text(json.dumps(stats, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(pf[pf.force_pn.astype(str) == ">=18"].T.to_string(), flush=True)
    print(ctrl.to_string(index=False), flush=True)
    print(pb.to_string(index=False), flush=True)
    print(en.to_string(index=False), flush=True)
    print(dsum.to_string(index=False), flush=True)
    print(amp.to_string(index=False), flush=True)
    print(pos.to_string(index=False), flush=True)
    print(json.dumps(coin), flush=True)

OVERLAP_SET: set = set()


def control_events(steps) -> pd.DataFrame:
    sm = steps[steps.cls.str.startswith("small") & steps.use].copy()
    lab = np.full(len(sm), "", dtype=object)
    lab[(sm.force_pn == 10.0).to_numpy()] = "10 pN (closed fork)"
    hi = (sm.force_pn >= 18).to_numpy()
    lab[hi & (sm.from_bp < OPEN_LO).to_numpy()] = ">=18 pN level < 50 bp"
    lab[hi & (sm.from_bp >= OPEN_HI).to_numpy()] = ">=18 pN level >= 900 bp"
    sm["control"] = lab
    return sm[sm.control != ""]


def control_comparison(steps, pls, injr, rng) -> pd.DataFrame:
    """Small steps per second, amplitude-matched in nm (>= 30 nm, < 100 bp B) and
    corrected by the recall of abrupt steps injected into the REAL traces of the
    same force: open fork at >= 18 pN vs the closed fork at 10 pN (no fork can
    move: below the ~15 pN unzipping force) and vs levels >= 900 bp B."""
    ev = injr.drop_duplicates(["i", "shape", "amp_nm", "t_inj_s"])
    ev = ev[ev["shape"] == "abrupt"]
    rec = {f: g.groupby("amp_nm").detected.mean() for f, g in ev.assign(F=ev.force_pn).groupby("F")}

    def rcl(f, a_nm):
        r = rec.get(f)
        if r is None:
            return np.full(np.size(a_nm), np.nan)
        return np.atleast_1d(np.interp(np.abs(a_nm), r.index.to_numpy(), r.to_numpy()))
    rows = []
    sm = steps[steps.cls.str.startswith("small") & steps.use]
    for s_ in SETS:
        P_ = pls[pls.set == s_]
        for lab, pp, ee, col in (
                ("10 pN, closed fork", P_[P_.force_pn == 10.0], sm[(sm.set == s_) & (sm.force_pn == 10.0)], "exp_total_s"),
                (">=18 pN, open fork", P_[(P_.force_pn >= 18) & P_.opened],
                 sm[(sm.set == s_) & (sm.force_pn >= 18) & sm.open_fork.astype(bool)], "exp_open_s"),
                (">=18 pN, level >= 900 bp B", P_[P_.force_pn >= 18],
                 sm[(sm.set == s_) & (sm.force_pn >= 18) & (sm.from_bp >= OPEN_HI)], "exp_full_s")):
            for amin in (0.0, 30.0, 40.0):
                e = ee[ee.step_nm.abs() >= amin]
                w = np.array([1.0 / max(rcl(f, a)[0], 0.1) for f, a in zip(e.force_pn, e.step_nm)]) if len(e) else np.array([])
                T = float(pp[col].sum())
                e2 = e.assign(w=w)
                pp2 = pp.assign(T=pp[col])
                num = e2.groupby("i").w.sum().reindex(pp2.i, fill_value=0).to_numpy()
                r, lo, hi_ = boot_ratio(num, pp2["T"].to_numpy(), tether_key(pp2), rng)
                rows.append({"set": s_, "population": lab, "min_abs_step_nm": amin, "time_s": round(T, 1),
                             "n_small": len(e), "n_fwd": int((e.cls == "small_fwd").sum()),
                             "n_back": int((e.cls == "small_back").sum()),
                             "rate_raw_per_s": round(len(e) / T, 4) if T > 0 else np.nan,
                             "rate_recall_corr_per_s": round(r, 4), "ci": f"{lo:.4f}-{hi_:.4f}",
                             "median_abs_step_nm": round(float(e.step_nm.abs().median()), 1) if len(e) else np.nan})
    return pd.DataFrame(rows)



def in_win(ev_t, lo, hi):
    return (ev_t > lo) & (ev_t <= hi)


def per_big(steps, pls, nev, rng) -> pd.DataFrame:
    """Per reference big event: P(>= 1 open-fork small jump within 1 s before / after),
    complete windows only; N2 null P on the same windows; far-baseline expectation."""
    rows = []
    hi_ = pls.set_index("i")
    tend = {i: plats_end for i, plats_end in zip(pls.i, pls.exp_total_s + T_SKIP)}
    sm = steps[steps.cls.str.startswith("small") & steps.use & steps.open_fork.astype(bool)]
    smg = {i: g.t_s.to_numpy() for i, g in sm.groupby("i")}
    nv = nev[(nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI) & (nev.null == "N2")]
    nvg = {i: g for i, g in nv.groupby("i")}
    big = steps[(steps.cls == "big_fwd") & (steps.force_pn >= 18)]
    events = [(r.i, r.t_s, "big jump", r.set, r.force_pn, r.tether, r.to_bp) for r in big.itertuples()]
    arr = pls[(pls.force_pn >= 18) & (pls.first_level_bp >= OPEN_LO) & (pls.first_level_bp < OPEN_HI)]
    events += [(r.i, 0.0, "arrival", r.set, r.force_pn, r.tether, r.first_level_bp) for r in arr.itertuples()]
    ev = []
    for i, tb, kind, s_, f, teth, land in events:
        te = tend[i]
        x = smg.get(i, np.array([]))
        g = nvg.get(i)
        rec = {"i": i, "kind": kind, "set": s_, "force_pn": f, "tether": teth, "landed_bp": land}
        for side, lo, hi in (("before", tb - WIN_S, tb), ("after", max(tb, T_SKIP) if kind == "arrival" else tb, tb + WIN_S)):
            if kind == "arrival" and side == "before":
                continue
            complete = lo >= T_SKIP - 1e-9 and hi <= te
            n_ = int(np.sum(in_win(x, lo, hi) if side == "after" else (x >= lo) & (x < hi)))
            pn = np.nan
            if g is not None:
                pn = float(np.mean([np.any(in_win(g[g.rep == r_].t_s.to_numpy(), lo, hi)) for r_ in range(N_NULL)]))
            else:
                pn = 0.0
            rec.update({f"{side}_complete": complete, f"{side}_n": n_, f"{side}_null_p": pn, f"{side}_len_s": hi - lo})
        ev.append(rec)
    ev = pd.DataFrame(ev)
    ev.to_csv(OUT / "per_big_event_windows.csv", index=False)
    # far-baseline rate
    for s_ in SETS:
        for kind in ("arrival", "big jump"):
            for f in FORCES + (">=18",):
                e = group_force(ev[(ev.set == s_) & (ev.kind == kind)], f)
                if not len(e):
                    continue
                P_ = group_force(pls[(pls.set == s_) & (pls.force_pn >= 18) & pls.opened], f)
                far_n = sm[sm.i.isin(P_.i) & (sm.wclass == "far")]
                far_T = float(P_["exp_open_far_s"].sum())
                rate_far = len(far_n) / far_T if far_T > 0 else np.nan
                rec = {"set": s_, "kind": kind, "force_pn": f, "n_events": len(e),
                       "rate_far_per_s": round(rate_far, 4) if np.isfinite(rate_far) else np.nan}
                for side in ("before", "after"):
                    if f"{side}_complete" not in e or (kind == "arrival" and side == "before"):
                        continue
                    c = e[e[f"{side}_complete"].astype(bool)]
                    if not len(c):
                        continue
                    pobs, lo_, hi_c = boot_ratio(c[f"{side}_n"] > 0, np.ones(len(c)), c.tether, rng)
                    pn = float(c[f"{side}_null_p"].mean())
                    L = float(c[f"{side}_len_s"].mean())
                    pfar = 1 - np.exp(-rate_far * L) if np.isfinite(rate_far) else np.nan
                    rec.update({f"n_{side}_complete": len(c), f"P_{side}": round(pobs, 3),
                                f"P_{side}_ci": f"{lo_:.3f}-{hi_c:.3f}", f"P_{side}_N2": round(pn, 3),
                                f"P_{side}_far_expect": round(pfar, 3), f"window_{side}_s": round(L, 2)})
                rows.append(rec)
    return pd.DataFrame(rows)


def enrichment(smo, pl18, nev, inj, rng) -> pd.DataFrame:
    """Rate of open-fork small jumps per s of open-fork time by window class
    (after arrival / after big / before big / 1-2 s / far), raw, N2-subtracted,
    recall-corrected (recall of injected 20-60 bp steps in the same window class)."""
    rows = []
    nv = nev[(nev.null == "N2") & (nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI)].copy()
    nv["wclass"] = window_class(nv.dt_nearest_big_s.to_numpy(), (nv.nearest_big_kind == "big jump").astype(int).to_numpy())
    inj = inj.copy()
    inj["wclass"] = window_class(inj.dt_nearest_big_s.to_numpy(), (inj.nearest_big_kind == "big jump").astype(int).to_numpy())
    for s_ in SETS:
        P_ = pl18[(pl18.set == s_) & pl18.opened]
        E_ = smo[smo.i.isin(P_.i)]
        N_ = nv[nv.i.isin(P_.i)]
        I_ = inj[(inj.set == s_) & inj.amp_bp.between(20, 60)]
        teth = P_.tether.to_numpy()
        far = None
        for w in WCLASSES:
            T = P_[f"exp_open_{w}_s"].to_numpy()
            n = E_[E_.wclass == w].groupby("i").size().reindex(P_.i, fill_value=0).to_numpy()
            r, lo, hi = boot_ratio(n, T, teth, rng)
            rn = len(N_[N_.wclass == w]) / N_NULL / T.sum() if T.sum() > 0 else np.nan
            rec = float(I_[I_.wclass == w].detected.mean()) if (I_.wclass == w).any() else np.nan
            row = {"set": s_, "window": w, "open_fork_s": round(float(T.sum()), 1), "n_small": int(n.sum()),
                   "rate_per_s": round(r, 4), "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                   "rate_N2_per_s": round(rn, 4), "recall_20_60bp": round(rec, 3),
                   "rate_net_recall_corr_per_s": round((r - rn) / rec, 4) if rec and np.isfinite(rec) else np.nan}
            rows.append(row)
        d = pd.DataFrame([x for x in rows if x["set"] == s_]).set_index("window")
        for w in WCLASSES:
            k = rows.index(next(x for x in rows if x["set"] == s_ and x["window"] == w))
            rows[k]["ratio_to_far_net_corr"] = round(d.loc[w, "rate_net_recall_corr_per_s"] /
                                                     d.loc["far", "rate_net_recall_corr_per_s"], 2) \
                if d.loc["far", "rate_net_recall_corr_per_s"] > 0 else np.nan
            # permutation-free Poisson test vs the far rate (raw counts)
            from scipy.stats import poisson
            lam_ = d.loc["far", "rate_per_s"] * d.loc[w, "open_fork_s"]
            rows[k]["p_poisson_gt_far"] = round(float(poisson.sf(d.loc[w, "n_small"] - 1, lam_)), 4) if lam_ > 0 else np.nan
    return pd.DataFrame(rows)


def dt_table(smo, pl18, nev, expo) -> pd.DataFrame:
    idx = {int(i): k for k, i in enumerate(expo["i"])}
    rows = []
    centers = 0.5 * (DT_EDGES[:-1] + DT_EDGES[1:])
    nv = nev[(nev.null == "N2") & (nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI)]
    for s_ in SETS:
        P_ = pl18[(pl18.set == s_) & pl18.opened]
        ks = [idx[int(i)] for i in P_.i if int(i) in idx]
        for kind, key in (("arrival", "dt_arrival"), ("big jump", "dt_big")):
            T = expo[key][ks].sum(0)
            e = smo[smo.i.isin(P_.i) & (smo.nearest_big_kind == kind)]
            h, _ = np.histogram(e.dt_nearest_big_s, DT_EDGES)
            nn = nv[nv.i.isin(P_.i) & (nv.nearest_big_kind == kind)]
            hn, _ = np.histogram(nn.dt_nearest_big_s, DT_EDGES)
            for c, n_, t_, m_ in zip(centers, h, T, hn):
                rows.append({"set": s_, "kind": kind, "dt_mid_s": round(c, 3), "n_small": int(n_),
                             "open_fork_s": round(float(t_), 2), "rate_per_s": n_ / t_ if t_ > 0 else np.nan,
                             "rate_N2_per_s": m_ / N_NULL / t_ if t_ > 0 else np.nan})
    return pd.DataFrame(rows)



AMP_BINS_DISC = (0.0, 25.0, 40.0, 60.0, 100.0)
INJ_GRID = np.array(INJ_AMPS_BP[:-1])          # 10-80 bp (small)


def discreteness(smo, nev, inj, injr, ctrl_ev, rng):
    """dBIC = BIC(continuous PL, 0-3 knots) - BIC(flat step) per small jump.
    Threshold theta per |step| bin = 95th percentile on N3 (continuous PL
    surrogate, ramps >= 0.1 s): 5 % false 'abrupt' by construction.
    Calibration of the same statistic:
      N3 / N2          surrogates (bootstrap noise);
      H1 boot          abrupt steps injected into the N2 surrogate;
      H1 real          abrupt steps injected into the REAL trace (real drift and
                       low-frequency noise) - the power used for the estimate;
      ramp 0.3 / 1.0 s linear ramps injected into the real trace (creep).
    Estimated abrupt fraction = (f_real - a_N3) / (power_H1real - a_N3)."""
    op = (nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI) & np.isfinite(nev.dbic_pl_minus_step)
    n3 = nev[op & (nev.null == "N3")].copy()
    n2 = nev[op & (nev.null == "N2")].copy()
    h1 = inj[inj.detected & inj.est_cls.astype(str).str.startswith("small") &
             np.isfinite(inj.dbic_pl_minus_step)].copy()
    h1["step_bp"] = h1.est_bp
    ir = injr[(injr.force_pn >= 18) & injr.detected & np.isfinite(injr.dbic_pl_minus_step)].copy()
    ir["step_bp"] = ir.piece_bp
    real = smo[np.isfinite(smo.dbic_pl_minus_step)].copy()
    ctrl = ctrl_ev[np.isfinite(ctrl_ev.dbic_pl_minus_step)].copy()
    for d in (n3, n2, real, h1, ir, ctrl):
        d["abin"] = np.digitize(d.step_bp.abs(), AMP_BINS_DISC[1:-1])
    th = {b: float(np.percentile(n3[n3.abin == b].dbic_pl_minus_step, 95)) if (n3.abin == b).sum() >= 20
          else float(np.percentile(n3.dbic_pl_minus_step, 95)) for b in range(len(AMP_BINS_DISC) - 1)}
    for d in (n3, n2, real, h1, ir, ctrl):
        d["theta"] = d.abin.map(th)
        d["discrete"] = d.dbic_pl_minus_step > d.theta
    near = lambda d: d.dt_nearest_big_s.abs() <= WIN_S
    ir["near"] = near(ir)
    real["near"] = near(real)
    hr = ir[ir["shape"] == "abrupt"]
    # matched on the DETECTED step size (|piece| bin; detected sizes near the
    # threshold are biased high), the force and near (<= 1 s) / far from big events
    hb = {k: g.discrete.to_numpy() for k, g in hr.groupby(["force_pn", "abin", "near"])}
    hbb = {k: g.discrete.to_numpy() for k, g in h1.groupby(["force_pn", "abin"])}
    pw_r, pw_b, h1s = [], [], []
    for r in real.itertuples():
        x = hb.get((r.force_pn, r.abin, bool(r.near)))
        if x is None or x.size < 8:
            x = hr[(hr.abin == r.abin) & (hr.near == bool(r.near))].discrete.to_numpy()
        pw_r.append(float(np.mean(x)) if x.size else np.nan)
        x2 = hbb.get((r.force_pn, r.abin))
        if x2 is None or x2.size < 8:
            x2 = h1[h1.abin == r.abin].discrete.to_numpy()
        pw_b.append(float(np.mean(x2)) if x2.size else np.nan)
        xs = hr[(hr.abin == r.abin) & np.isclose(hr.force_pn, r.force_pn)].dbic_pl_minus_step.to_numpy()
        h1s.append(float(rng.choice(xs)) if xs.size else np.nan)
    real["power_h1_real"] = pw_r
    real["power_h1_boot"] = pw_b
    real["h1_real_matched_dbic"] = h1s
    rows = []
    subsets = (("all", lambda d: np.ones(len(d), bool)),
               ("forward", lambda d: d.cls.astype(str).str.endswith("fwd")),
               ("backward", lambda d: d.cls.astype(str).str.endswith("back")),
               ("<= 1 s from a big event", lambda d: d.dt_nearest_big_s.abs() <= WIN_S),
               ("> 2 s from big events", lambda d: d.dt_nearest_big_s.abs() > FAR_S),
               ("|step| < 40 bp", lambda d: d.step_bp.abs() < 40),
               ("|step| 40-100 bp", lambda d: d.step_bp.abs() >= 40))
    for s_ in SETS:
        for lab, sel in subsets:
            r_ = real[(real.set == s_) & sel(real)]
            if len(r_) < 3:
                continue
            fo, lo, hi = boot_ratio(r_.discrete.astype(float), np.ones(len(r_)), tether_key(r_), rng)
            sub3 = n3[(n3.set == s_) & sel(n3)]
            a3 = float(sub3.discrete.mean()) if len(sub3) else 0.05
            irs = ir[ir.set == s_]
            # ramps: weighted to the real |step| bins of this subset
            wb = r_.abin.value_counts(normalize=True)
            ramp = {}
            for sh in ("ramp 0.3 s", "ramp 1.0 s"):
                rr_ = irs[irs["shape"] == sh].groupby("abin").discrete.mean()
                ramp[sh] = float(sum(wb.get(b, 0) * rr_.get(b, np.nan) for b in wb.index if b in rr_.index) /
                                 max(sum(wb.get(b, 0) for b in wb.index if b in rr_.index), 1e-9))
            pwr = float(np.nanmean(r_.power_h1_real))
            pi = (fo - a3) / (pwr - a3) if pwr > a3 else np.nan
            rows.append({"set": s_, "subset": lab, "n_real": len(r_),
                         "median_dBIC_real": round(float(r_.dbic_pl_minus_step.median()), 1),
                         "median_dBIC_N3": round(float(sub3.dbic_pl_minus_step.median()), 1) if len(sub3) else np.nan,
                         "median_dBIC_H1real_matched": round(float(np.nanmedian(r_.h1_real_matched_dbic)), 1),
                         "frac_real_abrupt": round(fo, 3), "ci": f"{lo:.3f}-{hi:.3f}",
                         "false_abrupt_N3": round(a3, 3),
                         "false_abrupt_ramp0.3s_real": round(float(ramp["ramp 0.3 s"]), 3),
                         "false_abrupt_ramp1.0s_real": round(float(ramp["ramp 1.0 s"]), 3),
                         "power_H1_boot": round(float(np.nanmean(r_.power_h1_boot)), 3),
                         "power_H1_real": round(pwr, 3),
                         "est_fraction_abrupt": round(float(np.clip(pi, 0, 1)), 2) if np.isfinite(pi) else np.nan})
    dsum = pd.DataFrame(rows)
    # calibration table: every population through the same statistic
    cal = []
    pops = [("real open fork, " + SET_LABEL[s_], real[real.set == s_]) for s_ in SETS]
    pops += [(f"control {c}", ctrl[ctrl.control == c]) for c in ctrl.control.unique()]
    pops += [("N2 noise surrogate", n2), ("N3 PL surrogate (ramps >= 0.1 s)", n3), ("H1 abrupt, bootstrap noise", h1)]
    pops += [(f"{sh}, real noise", ir[ir["shape"] == sh]) for sh in ("abrupt", "ramp 0.3 s", "ramp 1.0 s")]
    for lab, d in pops:
        if not len(d):
            continue
        cal.append({"population": lab, "n": len(d), "median_abs_step_bp": round(float(d.step_bp.abs().median()), 1),
                    "median_dBIC": round(float(d.dbic_pl_minus_step.median()), 1),
                    "frac_dBIC_gt_theta": round(float(d.discrete.mean()), 3),
                    "frac_dBIC_gt_6": round(float((d.dbic_pl_minus_step > 6).mean()), 3)})
    cal = pd.DataFrame(cal)
    cal.to_csv(OUT / "discreteness_calibration.csv", index=False)
    real.to_csv(OUT / "discreteness_real_events.csv", index=False)
    n3.to_csv(OUT / "discreteness_N3_events.csv.gz", index=False)
    ir.to_csv(OUT / "discreteness_injreal_pieces.csv.gz", index=False)
    print("theta (N3 95th pct by |step| bin):", {f"{AMP_BINS_DISC[k]:g}-{AMP_BINS_DISC[k + 1]:g}": round(v, 1)
                                                for k, v in th.items()}, flush=True)
    print(cal.to_string(index=False), flush=True)
    return real, dsum, cal, th


def _reflect_draw(rng, x, h, lo, hi, n):
    """Smooth null: Gaussian KDE (bandwidth h) of x, reflected into [lo, hi)."""
    v = rng.choice(x, n) + rng.normal(0.0, h, n)
    for _ in range(5):
        v = np.where(v < lo, 2 * lo - v, v)
        v = np.where(v >= hi, 2 * hi - v - 1e-9, v)
    return np.clip(v, lo, hi - 1e-9)


PERIODS = np.arange(8.0, 50.01, 0.5)
KDE_H = 15.0
AMP_LO = 15.0


def amplitude_stats(smo, nev, rng, n_null=2000):
    """Scale-free periodicity scan of small-jump amplitudes (bp B) in [15, 100):
    C(P) = <cos(2 pi a / P)>, P = 8-50 bp; null = smooth KDE (h = 15 bp, which
    removes >= 83 % of any modulation with P <= 50 bp) of the same amplitudes;
    z(P) against that null, global p from max z; power for quantised amplitudes
    (P0 = 10 / 20 / 40 bp, own measurement SE, +-10 % per-tether scale error)."""
    rows, scan = [], []
    for s_ in SETS:
        for lab, sel in (("forward", smo.cls == "small_fwd"), ("backward", smo.cls == "small_back"),
                         ("both |step|", np.ones(len(smo), bool))):
            d = smo[(smo.set == s_) & sel]
            a = d.step_bp.abs().to_numpy()
            se = d.se_bp.to_numpy()[a >= AMP_LO]
            a = a[a >= AMP_LO]
            if a.size < 10:
                continue
            n = a.size
            c_obs = np.array([np.mean(np.cos(2 * np.pi * a / P)) for P in PERIODS])
            nulls = np.array([[np.mean(np.cos(2 * np.pi * v / P)) for P in PERIODS]
                              for v in (_reflect_draw(rng, a, KDE_H, AMP_LO, BIG_BP, n) for _ in range(n_null))])
            mu, sd = nulls.mean(0), nulls.std(0)
            z = (c_obs - mu) / sd
            zmax_null = ((nulls - mu) / sd).max(1)
            j = int(np.argmax(z))
            p_glob = float((np.sum(zmax_null >= z[j]) + 1) / (zmax_null.size + 1))
            thr = float(np.percentile(zmax_null, 95))
            pw = {}
            for P0, scale_err, fq in ((10.0, 0.0, 1.0), (10.0, 0.1, 1.0), (20.0, 0.0, 1.0), (20.0, 0.1, 1.0),
                                      (40.0, 0.0, 1.0), (40.0, 0.1, 1.0), (20.0, 0.1, 0.5), (40.0, 0.1, 0.5)):
                if True:
                    hits = 0
                    for _ in range(300):
                        v = _reflect_draw(rng, a, KDE_H, AMP_LO, BIG_BP, n)
                        qm = rng.random(n) < fq
                        v = np.where(qm, np.maximum(np.round(v / P0), 1) * P0, v)
                        v = v * (1 + rng.normal(0, scale_err, n)) + rng.normal(0, 1, n) * rng.choice(se, n)
                        v = v[(v >= AMP_LO) & (v < BIG_BP)]
                        if v.size < 5:
                            continue
                        cz = (np.array([np.mean(np.cos(2 * np.pi * v / P)) for P in PERIODS]) - mu) / sd * np.sqrt(v.size / n)
                        hits += cz.max() > thr
                    pw[f"power_P{P0:.0f}{'_scale10pct' if scale_err else ''}{'_half_quantised' if fq < 1 else ''}"] = round(hits / 300, 3)
            # preferred sizes: 5 bp bins vs the smooth null
            edges = np.arange(AMP_LO, BIG_BP + 0.1, 5.0)
            h, _ = np.histogram(a, edges)
            hn = np.mean([np.histogram(_reflect_draw(rng, a, KDE_H, AMP_LO, BIG_BP, n), edges)[0] for _ in range(500)], 0)
            from scipy.stats import poisson
            pbin = poisson.sf(h - 1, np.maximum(hn, 1e-9))
            kbest = int(np.argmin(pbin))
            rows.append({"set": s_, "direction": lab, "n": n, "median_bp": round(float(np.median(a)), 1),
                         "iqr_bp": f"{np.percentile(a, 25):.1f}-{np.percentile(a, 75):.1f}",
                         "best_period_bp": float(PERIODS[j]), "z_max": round(float(z[j]), 2), "p_global": round(p_glob, 4),
                         "most_excess_bin_bp": f"{edges[kbest]:.0f}-{edges[kbest + 1]:.0f}",
                         "most_excess_obs_vs_null": f"{h[kbest]} vs {hn[kbest]:.1f}", "most_excess_p_local": round(float(pbin[kbest]), 4),
                         "median_se_bp": round(float(np.median(se)), 1), **pw})
            for P, zz, cc in zip(PERIODS, z, c_obs):
                scan.append({"set": s_, "direction": lab, "period_bp": P, "C_obs": round(float(cc), 4), "z": round(float(zz), 3),
                             "z_global_95": round(thr, 3)})
    return pd.DataFrame(rows), pd.DataFrame(scan)


SITE_ZONE = (290.0, 680.0)     # the 601 region +-25 bp: landmarks every 27-40 bp
NEAR_LM_BP = 10.0


def position_stats(smo, pl18, expo, inj, rng, n_mc=4000):
    """Start / end level (Luz_B) of open-fork small jumps vs a detector-aware
    occupancy null: open-fork dwell time per 10 bp bin of the same plateaus, each
    plateau weighted by the recall of a 30-40 bp step at its force and noise
    (injections). Landmark test inside the 601 zone (290-680 bp B), where the
    landmarks (entry, dyad-40, dyad, dyad+40, exit) are 27-40 bp apart: fraction
    within +-10 bp of a landmark and median distance, Monte Carlo vs the null
    restricted to the same zone."""
    from scipy.stats import chisquare, ks_2samp
    idx = {int(i): k for k, i in enumerate(expo["i"])}
    centers = 0.5 * (LUZ_EDGES[:-1] + LUZ_EDGES[1:])
    sb_edges = np.unique(np.percentile(inj.sigma_nm, [0, 25, 50, 75, 100]))
    inj = inj.assign(sbin=np.clip(np.searchsorted(sb_edges, inj.sigma_nm, side="right") - 1, 0, len(sb_edges) - 2))
    rec = inj[inj.amp_bp.isin([30.0, 40.0])].groupby(["force_pn", "sbin"]).detected.mean()
    rows, nul = [], []
    for s_ in SETS:
        P_ = pl18[(pl18.set == s_) & pl18.opened]
        sbin = np.clip(np.searchsorted(sb_edges, P_.sigma_nm, side="right") - 1, 0, len(sb_edges) - 2)
        w = np.array([rec.get((f, b), np.nan) for f, b in zip(P_.force_pn, sbin)])
        w = np.where(np.isfinite(w), w, np.nanmean(w))
        occ = np.zeros(centers.size)
        for i, wi in zip(P_.i, w):
            if int(i) in idx:
                occ += wi * expo["occupancy"][idx[int(i)]]
        dens = occ / occ.sum()
        for c, o in zip(centers, dens):
            nul.append({"set": s_, "luz_bp_B": c, "null_density": o})
        e = smo[smo.i.isin(P_.i)]
        zc = (centers >= SITE_ZONE[0]) & (centers < SITE_ZONE[1])
        dz = np.where(zc, dens, 0.0)
        dz = dz / dz.sum()
        for col in ("from_bp", "to_bp"):
            x = e[col].to_numpy()
            x = x[np.isfinite(x)]
            if x.size < 5:
                continue
            draws = rng.choice(centers, (200, x.size), p=dens) + rng.uniform(-5, 5, (200, x.size))
            reg_obs = pd.Series([SB.turn_region(v) for v in x]).value_counts()
            reg_null = pd.Series([SB.turn_region(v) for v in draws.ravel()]).value_counts(normalize=True)
            regs = sorted(set(reg_null.index) | set(reg_obs.index))
            o_ = np.array([reg_obs.get(r_, 0) for r_ in regs], float)
            e_ = np.array([reg_null.get(r_, 0) for r_ in regs], float) * x.size
            keep = e_ > 0
            chi = chisquare(o_[keep], e_[keep] * o_[keep].sum() / e_[keep].sum())
            ks = ks_2samp(x, draws.ravel())
            xz = x[(x >= SITE_ZONE[0]) & (x < SITE_ZONE[1])]
            dmin = lambda v: np.min(np.abs(v[..., None] - LM), axis=-1)
            if xz.size >= 5:
                dr = rng.choice(centers, (n_mc, xz.size), p=dz) + rng.uniform(-5, 5, (n_mc, xz.size))
                near_null = (dmin(dr) <= NEAR_LM_BP).mean(1)
                med_null = np.median(dmin(dr), 1)
                near_obs = float(np.mean(dmin(xz) <= NEAR_LM_BP))
                med_obs = float(np.median(dmin(xz)))
                p_near = float((np.sum(near_null >= near_obs) + 1) / (n_mc + 1))
                p_med = float((np.sum(med_null <= med_obs) + 1) / (n_mc + 1))
            else:
                near_obs = med_obs = p_near = p_med = np.nan
                near_null = med_null = np.array([np.nan])
            rows.append({"set": s_, "level": col, "n": x.size, "n_in_601_zone": int(xz.size),
                         "frac_within_10bp_of_landmark_obs": round(near_obs, 3),
                         "frac_within_10bp_null": round(float(np.nanmean(near_null)), 3), "p_enriched": round(p_near, 4),
                         "median_landmark_dist_obs_bp": round(med_obs, 1),
                         "median_landmark_dist_null_bp": round(float(np.nanmean(med_null)), 1), "p_closer": round(p_med, 4),
                         "ks_vs_occupancy_null_p": round(float(ks.pvalue), 4), "regions_chi2_p": round(float(chi.pvalue), 4),
                         "regions_obs": "; ".join(f"{r_} {int(v)}" for r_, v in zip(regs, o_)),
                         "regions_null_expect": "; ".join(f"{r_} {v:.1f}" for r_, v in zip(regs, e_ * o_.sum() / max(e_.sum(), 1)))})
    return pd.DataFrame(rows), pd.DataFrame(nul)


def relative_landing(smo) -> pd.DataFrame:
    rows = []
    for s_ in SETS:
        for kind in ("arrival", "big jump", "any"):
            d = smo[smo.set == s_]
            # the reference is the PREVIOUS big event (the level it landed on)
            d = d[np.isfinite(d.from_rel_landed_bp)]
            if kind != "any":
                d = d[(d.dt_since_prev_big_s.notna())]
                prev_is_arr = np.isclose(d.t_s - d.dt_since_prev_big_s, 0.0)
                d = d[prev_is_arr] if kind == "arrival" else d[~prev_is_arr]
            if not len(d):
                continue
            first = np.isclose(d.from_rel_landed_bp, 0.0, atol=1e-6)
            rows.append({"set": s_, "previous_big_event": kind, "n": len(d),
                         "frac_first_change_after_landing": round(float(first.mean()), 3),
                         "median_from_rel_bp": round(float(d.from_rel_landed_bp.median()), 1),
                         "median_to_rel_bp": round(float(d.to_rel_landed_bp.median()), 1),
                         "to_rel_iqr_bp": f"{d.to_rel_landed_bp.quantile(.25):.0f}-{d.to_rel_landed_bp.quantile(.75):.0f}",
                         "frac_to_rel_within_0_100": round(float(d.to_rel_landed_bp.between(0, 100).mean()), 3),
                         "median_dt_since_landing_s": round(float(d.dt_since_prev_big_s.median()), 2)})
    return pd.DataFrame(rows)



# ================================================================= figures (WHZ style)
C_BIG, C_FWD, C_BACK, C_RAW = "#E64B35", "#3C5488", "#00A087", "#8C8C8C"


def _export(fig, name):
    from whz_plot import export_panel
    export_panel(fig, OUT / name, formats=("pdf", "png"), dpi=300)


def _map601(ax, horizontal=True, lo=-1e9, hi=1e9):
    span = ax.axhspan if horizontal else ax.axvspan
    line = ax.axhline if horizontal else ax.axvline
    for a, b in (N1, N2):
        if b > lo and a < hi:
            span(a, b, color="#D9D9D9", alpha=0.45, lw=0, zorder=0)
    for dy in DYADS:
        for x, ls in ((dy, "--"), (dy - 40, ":"), (dy + 40, ":")):
            if lo <= x <= hi:
                line(x, color="#9E9E9E", lw=0.5, ls=ls, zorder=0)


def pick_examples(smo, pls, n_total=12):
    """User's 96 only, >= 18 pN, clean plateaus (sigma <= 30 nm): D026 R039 at 25.9 and
    26.9 pN (the user's screenshot), then the clearest small jumps (|step| / SE) after
    a big jump, before a big jump, after arrival and far from big events, plus the
    most abrupt-looking and the most ramp-like one; one panel per tether."""
    e = smo[(smo.set == "user96") & (smo.sigma_nm <= 30)].copy()
    e["z"] = e.step_bp.abs() / e.se_bp
    picks, seen = [], set()
    for f_ in (25.9, 26.9):
        r = e[(e.data_id == 26) & (e.bead == "R039") & np.isclose(e.force_pn, f_)].sort_values("t_s")
        if len(r):
            picks.append(("screenshot bead", r.iloc[-1]))
    seen.add((26, "R039"))
    for w, k in (("after big", 2), ("before big", 2), ("after arrival", 2), ("far", 2)):
        c = e[e.wclass == w].sort_values("z", ascending=False)
        n = 0
        for _, r in c.iterrows():
            if (r.data_id, r.bead) in seen or (w == "far" and np.isfinite(r.dbic_pl_minus_step) is False):
                continue
            picks.append((w, r))
            seen.add((r.data_id, r.bead))
            n += 1
            if n >= k:
                break
    for lab, c in (("most abrupt (dBIC)", e.sort_values("dbic_pl_minus_step", ascending=False)),
                   ("most ramp-like (dBIC)", e[e.step_bp.abs() >= 30].sort_values("dbic_pl_minus_step"))):
        for _, r in c.iterrows():
            if (r.data_id, r.bead) in seen:
                continue
            picks.append((lab, r))
            seen.add((r.data_id, r.bead))
            break
    return picks[:n_total]


def fig_examples(plats, smo, pls, steps):
    from whz_plot import raw_trace, fit_line, style_axes
    picks = pick_examples(smo, pls)
    nrow = int(np.ceil(len(picks) / 3))
    fig, axes = plt.subplots(nrow, 3, figsize=(7.09, 1.75 * nrow + 0.2), squeeze=False)
    rows = []
    for ax, (why, r) in zip(axes.flat, picks):
        q = plats[int(r.i)]
        f = q["force_pn"]
        t, y = q["t"], q["y"]
        ts, ys = search(q)
        sig = SB.sigma_diff(ys)
        st = fit_levels(ts, ys, sig, LAM)
        lz = luz_b(y, f)
        sp = steps[steps.i == r.i]
        t_c = float(r.t_s)
        bigs = sp[sp.cls == "big_fwd"].t_s.to_numpy()
        ref = np.concatenate([[0.0], bigs])
        tb = ref[np.argmin(np.abs(ref - t_c))]
        lo = max(0.0, min(t_c, tb) - 1.5)
        hi = min(t[-1], max(t_c, tb) + 1.5)
        if hi - lo < 4.0:
            lo, hi = max(0.0, 0.5 * (lo + hi) - 2.0), min(t[-1], 0.5 * (lo + hi) + 2.0)
        w = (t >= lo) & (t <= hi)
        raw_trace(ax, t[w], lz[w], color=C_RAW, filled=True, marker_area=2.0, edge_width=0.0, alpha=0.75,
                  rasterized=True)
        xs, yl = [], []
        for lv in st.itertuples():
            a_, b_ = max(lv.t0, lo), min(lv.t1, hi)
            if b_ > a_:
                xs += [a_, b_, np.nan]
                yl += [float(luz_b(lv.level_nm, f)[0])] * 2 + [np.nan]
        fit_line(ax, xs, yl, color="black", linewidth=0.9)
        dbs = []
        for c in sp.itertuples():
            if not lo <= c.t_s <= hi:
                continue
            col = C_BIG if c.cls == "big_fwd" else C_FWD if c.cls == "small_fwd" else C_BACK if c.cls == "small_back" else "#7E6148"
            fit_line(ax, [c.t_s, c.t_s], [c.from_bp, c.to_bp], color=col, linewidth=2.0, zorder=5)
            if c.cls.startswith("small") and np.isfinite(c.dbic_pl_minus_step):
                dbs.append(f"{c.dbic_pl_minus_step:+.0f}")
            rows.append({"panel": len(set(x["panel"] for x in rows)) if rows else 0, "why": why, "data_id": c.data_id,
                         "bead": c.bead, "force_pn": f, "t_s": round(c.t_s, 3), "class": c.cls,
                         "step_bp_B": round(c.step_bp, 1), "step_nm": round(c.step_nm, 1),
                         "from_bp_B": round(c.from_bp, 1), "to_bp_B": round(c.to_bp, 1),
                         "dt_nearest_big_s": round(c.dt_nearest_big_s, 3), "dBIC": c.dbic_pl_minus_step})
        if lo <= T_SKIP:
            ax.axvspan(0, T_SKIP, color="#EFEFEF", lw=0, zorder=0)
            ax.plot([0.0], [np.nanpercentile(lz[w], 2)], marker="^", ms=4, color=C_BIG, ls="none", zorder=6)
        yl_ = np.nanpercentile(lz[w], [0.5, 99.5])
        _map601(ax, True, yl_[0] - 30, yl_[1] + 30)
        ax.set_xlim(lo, hi)
        ax.set_ylim(yl_[0] - 30, yl_[1] + 30)
        ax.set_title(f"D{q['data_id']:03d} {q['bead']} {f:g} pN; {why}; dBIC {' '.join(dbs) or 'n/a'}",
                     fontsize=5.6, loc="left", pad=2)
        style_axes(ax, xlabel="Time in plateau (s)", ylabel="Luz B (bp)")
        ax.tick_params(labelsize=6)
        ax.xaxis.label.set_size(6.5)
        ax.yaxis.label.set_size(6.5)
    for ax in list(axes.flat)[len(picks):]:
        ax.remove()
    fig.tight_layout(h_pad=0.5, w_pad=0.6)
    _export(fig, "fig1_examples")
    pd.DataFrame(rows).to_csv(OUT / "fig1_examples_events.csv", index=False)


def fig_detection():
    from whz_plot import style_axes, fit_line
    se = pd.read_csv(OUT / "settle_changepoint_density.csv")
    rt = pd.read_csv(OUT / "detection_recall_by_force.csv")
    lim = pd.read_csv(OUT / "detection_limits_by_force.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.2))
    ax = axes[0, 0]
    for g, c, mk in (("10 pN", "#8491B4", "s"), (">=18 pN", "black", "o")):
        d = se[se.group == g]
        ax.plot(d.t_mid_s, d.mean_ext_minus_1p5_3s_median_nm, marker=mk, ms=3, ls="none", color=c,
                label=f"{g} ({int(d.n_plateaus.iloc[0])} plateaus)")
    ax.axvline(T_SKIP, color=C_BIG, lw=0.8, ls="--")
    ax.axhline(0, color="#BFBFBF", lw=0.6, zorder=0)
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    style_axes(ax, xlabel="Time since arrival (s)", ylabel="Extension - 1.5-3 s median (nm)", xlim=(0, 3))
    ax.set_title("a  settling after the 2 mm/s move", fontsize=8, loc="left")
    ax = axes[0, 1]
    d = se[se.group == ">=18 pN"]
    ax.stairs(d.small_change_points_per_plateau_per_s.to_numpy()[:-1] if False else d.small_change_points_per_plateau_per_s.to_numpy(),
              SETTLE_BINS, color="black", lw=1.0, label="small steps, >= 18 pN")
    ax.stairs(d.change_points_per_plateau_per_s.to_numpy(), SETTLE_BINS, color="#BFBFBF", lw=1.0, label="all steps, >= 18 pN")
    ax.axvline(T_SKIP, color=C_BIG, lw=0.8, ls="--")
    ax.legend(frameon=False, fontsize=6.5)
    style_axes(ax, xlabel="Time since arrival (s)", ylabel="Steps per plateau per s", xlim=(0, 3))
    ax.set_title("b  step density, OP on the whole plateau", fontsize=8, loc="left")
    ax = axes[1, 0]
    for f_, c in ((18.0, "#E64B35"), (30.0, "#3C5488"), (39.9, "#00A087")):
        r = rt[(rt.set == "user96") & np.isclose(rt.force_pn, f_)].sort_values("amp_bp")
        fit_line(ax, r.amp_bp, r.recall, color=c, linewidth=1.0, marker="o", ms=3, label=f"{f_:g} pN")
    ax.axhline(0.8, color="#BFBFBF", lw=0.6, ls=":")
    ax.axhline(0.5, color="#BFBFBF", lw=0.6, ls=":")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    style_axes(ax, xlabel="Injected step (bp, model B)", ylabel="Recall", xlim=(0, 125), ylim=(0, 1.02))
    ax.set_title("c  recall, own clock and noise (your 96)", fontsize=8, loc="left")
    ax = axes[1, 1]
    for s_, mk in (("user96", "o"), ("other_clean", "s")):
        l_ = lim[lim.set == s_].sort_values("force_pn")
        l_ = l_[l_.force_pn < 45]
        ax.plot(l_.force_pn, l_.A80_bp_B, marker=mk, ms=3.5, ls="none", color=SET_COL[s_], label=f"A80, {SET_LABEL[s_]}")
        ax.plot(l_.force_pn, l_.A50_bp_B, marker=mk, ms=3.5, ls="none", mfc="none", color=SET_COL[s_], label=f"A50, {SET_LABEL[s_]}")
    ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="lower center")
    style_axes(ax, xlabel="Plateau force (pN)", ylabel="Step size at recall (bp, model B)", ylim=(0, 70))
    ax.set_title("d  detectable size: 80 % / 50 % recall", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig2_settling_and_detection")


def fig_probability():
    from whz_plot import style_axes, fit_line
    pf = pd.read_csv(OUT / "small_jump_probability_by_force.csv")
    cm = pd.read_csv(OUT / "controls_amplitude_matched.csv")
    pf = pf[pf.force_pn.astype(str) != ">=18"].copy()
    pf["F"] = pf.force_pn.astype(float)
    pf = pf[pf.F < 45]
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6))
    ax = axes[0]
    for s_, mk, dx in (("user96", "o", -0.25), ("other_clean", "s", 0.25)):
        q = pf[pf.set == s_].sort_values("F")
        ax.errorbar(q.F + dx, q.P_plateau_ge1_small, yerr=[q.P_plateau_ge1_small - q.P_ci_lo, q.P_ci_hi - q.P_plateau_ge1_small],
                    fmt=mk, ms=3.5, color=SET_COL[s_], lw=0.8, capsize=1.5, label=SET_LABEL[s_])
    q = pf[pf.set == "user96"].sort_values("F")
    fit_line(ax, q.F, q.P_N2_null, color="black", linestyle="--", linewidth=0.9, label="N2 noise null")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    style_axes(ax, xlabel="Plateau force (pN)", ylabel="P(plateau has >= 1 small jump)", ylim=(0, 0.8))
    ax.set_title("a  per opened plateau", fontsize=8, loc="left")
    ax = axes[1]
    for s_, mk, dx in (("user96", "o", -0.25), ("other_clean", "s", 0.25)):
        q = pf[pf.set == s_].sort_values("F")
        ax.errorbar(q.F + dx, q.rate_small_per_s, yerr=[q.rate_small_per_s - q.rate_ci_lo, q.rate_ci_hi - q.rate_small_per_s],
                    fmt=mk, ms=3.5, color=SET_COL[s_], lw=0.8, capsize=1.5, label=f"all, {SET_LABEL[s_]}")
    c10 = cm[(cm.set == "user96") & (cm.population == "10 pN, closed fork") & (cm.min_abs_step_nm == 0)].rate_raw_per_s.iloc[0]
    ax.axhline(c10, color=SET_COL["user96"], lw=0.8, ls=":", label="10 pN closed fork, your 96")
    ax.legend(frameon=False, fontsize=6.0, loc="upper right")
    style_axes(ax, xlabel="Plateau force (pN)", ylabel="Small jumps per s open fork", ylim=(0, 0.21))
    ax.set_title("b  all detected small jumps", fontsize=8, loc="left")
    ax = axes[2]
    q = pf[pf.set == "user96"].sort_values("F")
    ax.plot(q.F, q.rate_detectable_net_per_s, "o", ms=3.5, color=SET_COL["user96"], label=">= A80, recall-corr., your 96")
    q2 = pf[pf.set == "other_clean"].sort_values("F")
    ax.plot(q2.F, q2.rate_detectable_net_per_s, "s", ms=3.5, color=SET_COL["other_clean"], label=">= A80, other clean")
    c40 = cm[(cm.set == "user96") & (cm.population == "10 pN, closed fork") & (cm.min_abs_step_nm == 40)].rate_recall_corr_per_s.iloc[0]
    ax.axhline(c40, color=SET_COL["user96"], lw=0.8, ls=":", label="10 pN control, >= 40 nm")
    ax.legend(frameon=False, fontsize=6.0, loc="upper left")
    style_axes(ax, xlabel="Plateau force (pN)", ylabel="Small jumps per s open fork", ylim=(0, 0.07))
    ax.set_title("c  detectable range, 49-100 bp", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig3_probability_vs_force")


def fig_amplitudes():
    from whz_plot import style_axes, fit_line
    smo = pd.read_csv(OUT / "small_jumps_open_fork.csv")
    nev = pd.read_csv(OUT / "null_small_jumps.csv.gz")
    sc = pd.read_csv(OUT / "periodicity_scan.csv")
    st = pd.read_csv(OUT / "steps_all_flagged.csv.gz")
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6))
    ax = axes[0]
    edges = np.arange(0, 101, 5.0)
    u = smo[smo.set == "user96"]
    ax.stairs(np.histogram(u[u.cls == "small_fwd"].step_bp.abs(), edges)[0], edges, color=C_FWD, lw=1.1,
              label=f"forward (n = {int((u.cls == 'small_fwd').sum())})")
    ax.stairs(np.histogram(u[u.cls == "small_back"].step_bp.abs(), edges)[0], edges, color=C_BACK, lw=1.1,
              label=f"backward (n = {int((u.cls == 'small_back').sum())})")
    op = (nev.from_bp >= OPEN_LO) & (nev.from_bp < OPEN_HI) & (nev.set == "user96")
    h3 = np.histogram(nev[op & (nev.null == "N3")].step_bp.abs(), edges)[0] / N_NULL
    ax.stairs(h3, edges, color="black", ls="--", lw=0.9, label="N3 creep surrogate")
    for x in (10, 20, 40, 80):
        ax.axvline(x, color="#BFBFBF", lw=0.5, ls=":", zorder=0)
    ax.legend(frameon=False, fontsize=6.0)
    style_axes(ax, xlabel="|Small jump| (bp, model B)", ylabel="Small jumps per 5 bp", xlim=(0, 100))
    ax.set_title("a  amplitudes, your 96, open fork", fontsize=8, loc="left")
    ax = axes[1]
    en = np.arange(0, 131, 10.0)
    c10 = st[(st.set == "user96") & (st.force_pn == 10.0) & st.cls.str.startswith("small") & st.use]
    T10 = pd.read_csv(OUT / "plateau_summary.csv")
    T10 = T10[(T10.set == "user96") & (T10.force_pn == 10.0)].exp_total_s.sum()
    To = pd.read_csv(OUT / "plateau_summary.csv")
    To = To[(To.set == "user96") & (To.force_pn >= 18) & To.opened].exp_open_s.sum()
    ax.stairs(np.histogram(u.step_nm.abs(), en)[0] / To * 1000, en, color=C_FWD, lw=1.1, label="open fork, >= 18 pN")
    ax.stairs(np.histogram(c10.step_nm.abs(), en)[0] / T10 * 1000, en, color="#8491B4", lw=1.1, ls="--",
              label="closed fork, 10 pN")
    ax.legend(frameon=False, fontsize=6.0)
    style_axes(ax, xlabel="|Small jump| (nm)", ylabel="Per 1000 s per 10 nm", xlim=(0, 130))
    ax.set_title("b  vs no-fork control (raw)", fontsize=8, loc="left")
    ax = axes[2]
    for s_, d_, c, ls in (("user96", "forward", C_FWD, "-"), ("user96", "both |step|", "black", "-"),
                          ("other_clean", "forward", SET_COL["other_clean"], ":")):
        q = sc[(sc.set == s_) & (sc.direction == d_)]
        fit_line(ax, q.period_bp, q.z, color=c, linestyle=ls, linewidth=0.9, label=f"{SET_LABEL[s_]}, {d_}")
    thr = sc[(sc.set == "user96") & (sc.direction == "forward")].z_global_95.iloc[0]
    ax.axhline(thr, color=C_BIG, lw=0.8, ls="--", label="global 5 % level")
    ax.legend(frameon=False, fontsize=5.8, loc="lower right")
    style_axes(ax, xlabel="Period P (bp, model B)", ylabel="z of <cos(2 pi a / P)>", xlim=(8, 50), ylim=(-4, 5))
    ax.set_title("c  periodicity scan", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig4_amplitudes")


def fig_timing():
    from whz_plot import style_axes
    dt = pd.read_csv(OUT / "rate_vs_dt_to_big.csv")
    en = pd.read_csv(OUT / "enrichment_by_window.csv")
    pb = pd.read_csv(OUT / "per_big_event.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.6))
    far = float(en[(en.set == "user96") & (en.window == "far")].rate_per_s.iloc[0])
    for ax, kind, ttl in ((axes[0], "big jump", "a  around big jumps (>= 100 bp)"), (axes[1], "arrival", "b  after arrival (t = 0)")):
        d = dt[(dt.set == "user96") & (dt.kind == kind)].sort_values("dt_mid_s")
        # 0.2 s bins for readability
        e2 = np.round(np.arange(-3.0, 3.0001, 0.2), 3)
        n_ = np.histogram(d.dt_mid_s, e2, weights=d.n_small)[0]
        T_ = np.histogram(d.dt_mid_s, e2, weights=d.open_fork_s)[0]
        nn = np.histogram(d.dt_mid_s, e2, weights=d.rate_N2_per_s * d.open_fork_s)[0]
        with np.errstate(invalid="ignore", divide="ignore"):
            r_ = np.where(T_ > 10, n_ / T_, np.nan)
            rn = np.where(T_ > 10, nn / T_, np.nan)
        ok = np.isfinite(r_)
        x = 0.5 * (e2[:-1] + e2[1:])
        err = np.where(T_ > 10, np.sqrt(np.maximum(n_, 1)) / np.maximum(T_, 1e-9), np.nan)
        ax.errorbar(x[ok], r_[ok], yerr=err[ok], fmt="o", ms=3, color=C_FWD, lw=0.7, capsize=1.2, label="your 96")
        ax.stairs(np.nan_to_num(rn), e2, color="black", ls=":", lw=0.9, label="N2 noise null")
        ax.axhline(far, color="black", lw=0.8, ls="--", label="far (> 2 s) rate")
        ax.axvline(0, color=C_BIG, lw=0.8)
        ax.legend(frameon=False, fontsize=6.0, loc="upper right")
        style_axes(ax, xlabel=f"t - t({'big jump' if kind == 'big jump' else 'arrival'}) (s)",
                   ylabel="Small jumps per s open fork", xlim=(-3, 3) if kind == "big jump" else (0, 3), ylim=(0, 0.6))
        ax.set_title(ttl, fontsize=8, loc="left")
    ax = axes[2]
    q = pb[(pb.set == "user96") & (pb.force_pn.astype(str) == ">=18")]
    labels, obs, lo, hi, exp_, n2 = [], [], [], [], [], []
    for kind, side in (("arrival", "after"), ("big jump", "before"), ("big jump", "after")):
        r = q[q.kind == kind].iloc[0]
        labels.append(f"{side} {kind}")
        obs.append(r[f"P_{side}"])
        a_, b_ = [float(v) for v in str(r[f"P_{side}_ci"]).split("-")]
        lo.append(r[f"P_{side}"] - a_)
        hi.append(b_ - r[f"P_{side}"])
        exp_.append(r[f"P_{side}_far_expect"])
        n2.append(r[f"P_{side}_N2"])
    x = np.arange(3)
    ax.errorbar(x, obs, yerr=[lo, hi], fmt="o", ms=4, color=C_FWD, lw=0.8, capsize=2, label="observed")
    ax.plot(x, exp_, "_", ms=12, mew=1.5, color="black", label="far-rate expectation")
    ax.plot(x, n2, "x", ms=4, color="black", label="N2 noise null")
    ax.set_xticks(x)
    ax.set_xticklabels(["0.3-1 s after\narrival", "1 s before\nbig jump", "1 s after\nbig jump"], fontsize=6)
    ax.legend(frameon=False, fontsize=6.0, loc="upper left")
    style_axes(ax, xlabel="Window", ylabel="P(>= 1 small jump)", xlim=(-0.5, 2.5), ylim=(0, 0.45))
    ax.set_title("c  per reference event, your 96", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig5_timing_vs_big_jumps")


def fig_positions():
    from whz_plot import style_axes, fit_line
    smo = pd.read_csv(OUT / "small_jumps_open_fork.csv")
    nul = pd.read_csv(OUT / "position_occupancy_null.csv")
    fig, axes = plt.subplots(2, 1, figsize=(7.09, 4.6))
    edges = np.arange(-100.0, 1001.0, 20.0)
    for ax, s_ in ((axes[0], "user96"), (axes[1], "other_clean")):
        u = smo[smo.set == s_]
        _map601(ax, horizontal=False)
        h, _ = np.histogram(u.from_bp, edges)
        ax.stairs(h, edges, color=C_FWD, lw=1.1, label=f"start level (n = {len(u)})")
        h2, _ = np.histogram(u.to_bp, edges)
        ax.stairs(h2, edges, color=C_BIG, lw=0.9, ls="--", label="end level")
        nn = nul[nul.set == s_]
        dens = np.histogram(nn.luz_bp_B, edges, weights=nn.null_density)[0]
        fit_line(ax, 0.5 * (edges[:-1] + edges[1:]), dens * len(u), color="black", linewidth=0.9,
                 label="occupancy null (recall-weighted)")
        ax.axvline(971, color="black", lw=0.6, ls="-.")
        ax.legend(frameon=False, fontsize=6.5, loc="upper right")
        style_axes(ax, xlabel="Luz, model B (bp from the 10 pN zero)", ylabel="Small jumps per 20 bp", xlim=(-100, 1000))
        ax.set_title(f"{'a' if s_ == 'user96' else 'b'}  {SET_LABEL[s_]}: 601 sites grey, dyad dashed, dyad +-40 dotted",
                     fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig6_positions_601_map")
    # relative to the level the previous big event landed on
    fig, ax = plt.subplots(1, 1, figsize=(3.39, 2.65))
    u = smo[smo.set == "user96"]
    e = np.arange(-100, 301, 10.0)
    prev_arr = np.isclose(u.t_s - u.dt_since_prev_big_s, 0.0)
    ax.stairs(np.histogram(u[prev_arr].to_rel_landed_bp, e)[0], e, color=C_FWD, lw=1.1, label="after arrival")
    ax.stairs(np.histogram(u[~prev_arr].to_rel_landed_bp, e)[0], e, color=C_BIG, lw=1.1, label="after a big jump")
    ax.axvline(0, color="#BFBFBF", lw=0.6)
    for x in (40, 80):
        ax.axvline(x, color="#BFBFBF", lw=0.5, ls=":")
    ax.legend(frameon=False, fontsize=6.5)
    style_axes(ax, xlabel="End level - landing level (bp, model B)", ylabel="Small jumps per 10 bp", xlim=(-100, 300))
    fig.tight_layout()
    _export(fig, "fig6b_relative_to_landing")


def fig_discreteness():
    from whz_plot import style_axes, fit_line
    real = pd.read_csv(OUT / "discreteness_real_events.csv")
    n3 = pd.read_csv(OUT / "discreteness_N3_events.csv.gz")
    ir = pd.read_csv(OUT / "discreteness_injreal_pieces.csv.gz")
    st = pd.read_csv(OUT / "steps_all_flagged.csv.gz")
    cal = pd.read_csv(OUT / "discreteness_calibration.csv")
    c10 = st[(st.force_pn == 10.0) & st.cls.str.startswith("small") & st.use & np.isfinite(st.dbic_pl_minus_step)]
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.8), gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]

    def ecdf(ax, x, **kw):
        x = np.sort(np.asarray(x, float))
        x = x[np.isfinite(x)]
        fit_line(ax, x, np.arange(1, x.size + 1) / x.size, **kw)
    ecdf(ax, real[real.set == "user96"].dbic_pl_minus_step, color=C_FWD, linewidth=1.3, label="your 96, open fork")
    ecdf(ax, c10.dbic_pl_minus_step, color="#8491B4", linewidth=1.0, linestyle="--", label="10 pN closed fork")
    rng = np.random.default_rng(8)
    ru = real[real.set == "user96"]
    for sh, c, ls, lab_ in (("abrupt", C_BIG, "-", "abrupt step, size-matched"),
                            ("ramp 1.0 s", "black", ":", "1 s ramp, size-matched")):
        g = {k: v.dbic_pl_minus_step.to_numpy() for k, v in ir[ir["shape"] == sh].groupby(["force_pn", "abin"])}
        smp = []
        for r in ru.itertuples():
            x = g.get((r.force_pn, r.abin))
            if x is not None and x.size:
                smp.extend(rng.choice(x, 20))
        ecdf(ax, smp, color=c, linewidth=1.0, linestyle=ls, label=lab_)
    ax.axvline(0, color="#BFBFBF", lw=0.6)
    ax.legend(frameon=False, fontsize=6.0, loc="upper left")
    style_axes(ax, xlabel="dBIC = BIC(continuous) - BIC(step)", ylabel="Cumulative fraction", xlim=(-40, 25), ylim=(0, 1.01))
    ax.set_title("a  same statistic on every population", fontsize=8, loc="left")
    ax = axes[1]
    order = ["N3 PL surrogate (ramps >= 0.1 s)", "ramp 1.0 s, real noise", "ramp 0.3 s, real noise", "abrupt, real noise",
             "control 10 pN (closed fork)", "real open fork, your 96", "real open fork, other clean"]
    lab = ["N3 PL", "ramp 1 s", "ramp 0.3 s", "abrupt", "10 pN ctrl", "your 96", "other clean"]
    v = [float(cal[cal.population == o].frac_dBIC_gt_theta.iloc[0]) if (cal.population == o).any() else np.nan for o in order]
    cols = ["black", "black", "black", C_BIG, "#8491B4", C_FWD, SET_COL["other_clean"]]
    for k, (vv, c) in enumerate(zip(v, cols)):
        ax.plot([k], [vv], "o", ms=4.5, color=c)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=6)
    ax.axhline(0.05, color="#BFBFBF", lw=0.6, ls=":")
    style_axes(ax, xlabel="Population", ylabel="Fraction called abrupt", ylim=(0, 0.7), xlim=(-0.5, len(order) - 0.5))
    ax.set_title("b  dBIC > N3 95th percentile", fontsize=8, loc="left")
    fig.tight_layout()
    _export(fig, "fig7_discreteness")


def stage_figs(plats, qc, data):
    from whz_plot import apply_style
    apply_style()
    smo = pd.read_csv(OUT / "small_jumps_open_fork.csv")
    pls = pd.read_csv(OUT / "plateau_summary.csv")
    steps = pd.read_csv(OUT / "steps_all_flagged.csv.gz")
    for fn, a in ((fig_examples, (plats, smo, pls, steps)), (fig_detection, ()), (fig_probability, ()),
                  (fig_amplitudes, ()), (fig_timing, ()), (fig_positions, ()), (fig_discreteness, ())):
        fn(*a)
        print(f"  {fn.__name__} done", flush=True)


# ================================================================= driver (stages)
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--recollect", action="store_true")
    ap.add_argument("--stage", default="all")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    data = load_cache(args.recollect)
    plats = data["plateaus"]
    _init(plats)
    print(f"{len(plats)} plateaus loaded ({time.time() - t0:.0f} s)", flush=True)
    todo = STAGES if args.stage == "all" else (STAGES[STAGES.index(args.stage[:-1]):]
                                                if args.stage.endswith("+") else (args.stage,))
    qc = stage_qc(plats, data) if "qc" in todo else pd.read_csv(OUT / "plateaus_qc.csv")
    if "settle" in todo:
        stage_settle(plats, qc)
        print(f"[{time.time() - t0:.0f} s] settle", flush=True)
    if "detect" in todo:
        stage_detect(plats, qc)
        print(f"[{time.time() - t0:.0f} s] detect", flush=True)
    if "nulls" in todo:
        stage_nulls(plats, qc)
        print(f"[{time.time() - t0:.0f} s] nulls", flush=True)
    if "inject" in todo:
        stage_inject(plats, qc)
        print(f"[{time.time() - t0:.0f} s] inject", flush=True)
    if "injreal" in todo:
        stage_injreal(plats, qc)
        print(f"[{time.time() - t0:.0f} s] injreal", flush=True)
    if "stats" in todo:
        stage_stats(plats, qc, data)
        print(f"[{time.time() - t0:.0f} s] stats", flush=True)
    if "figs" in todo:
        stage_figs(plats, qc, data)
        print(f"[{time.time() - t0:.0f} s] figs", flush=True)
    return 0


STAGES = ("qc", "settle", "detect", "nulls", "inject", "injreal", "stats", "figs")

if __name__ == "__main__":
    raise SystemExit(main())
