#!/usr/bin/env python3
"""20261009 force-titration batch: the whole force ladder inside ONE file.

Same construct as 20261004 (YTrunk_2601_980bp, two Widom 601 sites at 314-460
and 511-657 bp, normalised to 980 bp) and the same magnet/bead setup, but a
different protocol. On 20261004 each file was one ramp to one Fmax and the
ladder was built across files of a field. Here every titration file runs a
G-code step script (the H5 `Steps`/`Script`):

    rest 5.00 mm (15 s) -> 2.50 -> 1.70 -> 1.20 -> 1.15 -> 1.10 -> 1.00
    -> 0.95 -> 0.90 -> 0.85 -> 0.80 mm, each held 10 s, back to 5.00 mm
    for 15 s between steps, magnet moving at 2 mm/s (F120 mm/min).

data_003 is the first version of the script: 15 s holds, 20 s rests and an
extra 0.70 mm step at the end. The short files (23-33 s, magnet 2.5-17.5 mm,
no step script) are low-force checks taken before a titration in the same
field; they are listed in the manifest and not analysed.

So one bead = one file = one ladder, and every step is a cycle of its own:
up move (5.00 -> position), 10 s plateau, down move, 15 s rest. The 20261004
rules carry over with the step taking the place of the file:

- Zero: from each step's own up move at 14 pN (SOP anchor, raw samples,
  pulled below any step in its window). The 2.50 mm step (10 pN) never
  reaches 14 pN; it takes the zero of the nearest step in time that has one,
  and `zero_source` says which.
- Length scale: ONE per bead, from its highest step at >= 38.6 pN (0.80 mm =
  39.9 pN; 0.70 mm = 46 pN in data_003) in which the bead is tracked to the
  end of the plateau. Full anchor = median extension over the last second of
  that plateau, at the plateau force: the step protocol's down move crosses
  38 +- 4 pN in 1-3 samples, too few for the SOP's down-branch estimator,
  while the plateau samples the same fully loaded state for 10 s. Bead kind from that span /
  797 nm (dWLC span of a fully open 980 bp construct): full 0.80-1.25,
  partial < 0.80, over-long > 1.25. Partial and over-long beads take the
  file's full-unzip median span (>= 3 full beads), else the batch median.
- Tracking: tracking_quality.file_tracking on every file (20260913 criteria);
  no bead / lost before the pull / LUT level jump > 1.5 um -> not normalised;
  lost -> samples after the loss masked.
- Stall per step = median Luz over the last second of the first 10 s of the
  plateau (= the last second of a 10 s hold; for data_003's 15 s holds the
  9-10 s window, so all files compare at the same time after the step; the
  end-of-hold value is kept as `luz_end_of_hold_bp`). Used only if the hold
  lasted >= 9 s, the bead is tracked to the end of that window, the fork is
  stationary in it (IQR <= 100 bp) and the value is within -150..1470 bp.

Forces use the calibration the user gave for 20261004 (1.7 / 1.0 / 0.7 mm =
18 / 30 / 46 pN, piecewise log-linear); positions outside 0.7-1.7 mm (2.50 mm
-> 10.0 pN, 5.00 mm rest -> 1.6 pN) are extrapolations of the end segments.
ASSUMPTION (2026-10-11): same magnets and beads as 20261004, so the same
calibration; stated in the README and the manifest.

    ~/.venvs/trace_editor_bm/bin/python process_20261009.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from whz_plot import PALETTE, apply_style, style_axes

from process_tierA_sop import drift_corrected_nm, force_from_shift
from reprocess_26_core_sop import ANCHOR_HALF_WINDOW, estimate_anchor, first_step_force
from whz_all_934 import (RAW_FULL_MAX_REACH_PN, RAW_FULL_MIN_POINTS,
                         RAW_ZERO_MAX_REACH_PN, RAW_ZERO_MIN_POINTS)
from tracking_quality import file_tracking
from process_20261004 import (
    CONSTRUCT_BP, F_OWN_STANDARD, FULL_UNZIP_RATIO, LEVEL_JUMP_UM, MAX_SPAN_NM,
    MIN_FULL_BEADS_PER_FOV, MIN_SPAN_NM, NUCLEOSOME_SITES_BP, expected_span_nm, normalise,
)
from tether_stall_titration_20261004 import (
    END_REACH_S, END_WINDOW_S, MAX_END_IQR_BP, MOVING_SLOPE_BP_S, REQUIRED_PLATEAU_S,
    STALL_GUARD_BP, frame_dt, level_colours,
)

RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/20261009")
DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/"
            "Proecessed Traces/WHZ_20261009")
CALIBRATION = "user_20261004"
REST_MM = 5.00
POS_TOL_MM = 0.004          # a sample is "at" a position within this
ZERO_TARGET_PN, FULL_TARGET_PN = 14.0, 38.0
STALL_AT_S = 10.0           # compare every file 10 s after the step
# Zero: one estimate per step from a 2 mm/s move scatters by ~40-90 nm (median
# 37 nm, 25 % > 91 nm step to step on this batch), because the move crosses
# 14 +- 2.5 pN in a handful of frames. The handle state does not change between
# steps, so each bead's zeros are pooled into a Theil-Sen line against time
# (it follows the residual drift); a step whose own zero sits > ZERO_OUTLIER_NM
# off that line was on another tracking level and is not used.
ZERO_OUTLIER_NM = 150.0
MIN_ZERO_STEPS = 3
FULL_MAX_IQR_NM = 100.0     # reference plateau must be stationary (no level flipping)
NEVER_UNZIPS_BP = 50.0      # 20261004 rule: reported, not paired
SITE1_ENTRY, SITE2_ENTRY = NUCLEOSOME_SITES_BP[0][0], NUCLEOSOME_SITES_BP[1][0]
FULL_OPEN_BP = 0.92 * CONSTRUCT_BP


# ---------------------------------------------------------------------------
# reading and segmenting
# ---------------------------------------------------------------------------
def titration_files() -> tuple[list[int], list[dict]]:
    """Files with a step script are titrations; the rest are listed, not used."""
    titr, others = [], []
    luts = {p.stem for p in RAW.glob("*.lut")}
    for p in sorted(RAW.glob("data_*.h5")):
        d = int(p.stem.split("_")[1])
        with h5py.File(p) as f:
            has_steps = "Steps" in f
            shift = (np.asarray(f["Traces/Shift (mm)"]) if "Shift (mm)" in f["Traces"]
                     else np.array([np.nan]))
            t = np.asarray(f["Traces/Time (s)"])
        if has_steps:
            titr.append(d)
        else:
            others.append({"data_id": d, "kind": "LUT" if p.stem in luts else "low-force check",
                           "duration_s": round(float(t[-1] - t[0]), 1),
                           "shift_mm": [round(float(np.nanmin(shift)), 2) if np.isfinite(shift).any() else None,
                                        round(float(np.nanmax(shift)), 2) if np.isfinite(shift).any() else None]})
    return titr, others


def load(d: int) -> dict:
    with h5py.File(RAW / f"data_{d:03d}.h5") as f:
        tr = f["Traces"]
        t = np.asarray(tr["Time (s)"], float)
        s = np.asarray(tr["Shift (mm)"], float)
        rois = sorted(tr["Z (um)"].keys())
        z = {r: np.asarray(tr["Z (um)"][r], float) for r in rois}
        steps = f["Steps"][()]
        xy = f["Parameters"][()]
    targets = [float(p) for p in steps["Position (mm)"] if abs(p - REST_MM) > 1e-6]
    return {"data_id": d, "time_s": t, "shift_mm": s, "force_pn": force_from_shift(s, CALIBRATION),
            "rois": rois, "z_um": z, "targets": targets,
            "xy": {f"{int(r['ROI']):03d}": (float(r["X (pix)"]), float(r["Y (pix)"])) for r in xy}}


def runs(mask: np.ndarray) -> list[np.ndarray]:
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    return np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)


def segment_steps(run: dict) -> list[dict]:
    """Walk the script in order: rest, up move, plateau, down move per step."""
    s, n = run["shift_mm"], run["shift_mm"].size
    rest = np.abs(s - REST_MM) <= POS_TOL_MM
    out, cursor = [], 0
    for k, pos in enumerate(run["targets"], start=1):
        at = np.abs(s - pos) <= POS_TOL_MM
        cand = [r for r in runs(at) if r[0] >= cursor]
        if not cand:
            break
        pl = max(cand[:3], key=len) if len(cand) > 1 else cand[0]   # the hold, not a pass-through
        before = np.flatnonzero(rest[:pl[0]])
        up_from = before[-1] + 1 if before.size else cursor
        after = np.flatnonzero(rest[pl[-1] + 1:])
        down_to = pl[-1] + 1 + after[0] if after.size else n
        out.append({"step": k, "pos_mm": pos,
                    "force_pn": float(force_from_shift(np.array([pos]), CALIBRATION)[0]),
                    "up": np.arange(up_from, pl[0]), "plateau": pl,
                    "down": np.arange(pl[-1] + 1, down_to)})
        cursor = down_to
    return out


# ---------------------------------------------------------------------------
# anchors and normalisation
# ---------------------------------------------------------------------------
def branch_frame(force, ext, branch):
    return pd.DataFrame({"force_bin_pn": force, "extension_raw_nm": ext, "branch": branch})


def step_anchors(run: dict, steps: list[dict], roi: str, ext: np.ndarray) -> list[dict]:
    t, f = run["time_s"], run["force_pn"]
    rows = []
    for st in steps:
        rec = {"step": st["step"], "pos_mm": st["pos_mm"], "plateau_force_pn": round(st["force_pn"], 2),
               "t_plateau_start_s": round(float(t[st["plateau"][0]] - t[0]), 2),
               "t_plateau_end_s": round(float(t[st["plateau"][-1]] - t[0]), 2),
               "hold_s": round(float(t[st["plateau"][-1]] - t[st["plateau"][0]]), 2),
               "z_zero_raw_nm": np.nan, "zero_force_pn": np.nan, "zero_note": "",
               "z_full_raw_nm": np.nan, "full_force_pn": np.nan, "t_full_done_s": np.nan}
        up = st["up"]
        if up.size and f[up].max() >= ZERO_TARGET_PN + ANCHOR_HALF_WINDOW:
            zbr = branch_frame(f[up], ext[up], "up")
            step_f = first_step_force(zbr, 10.0, ZERO_TARGET_PN, lag=5)
            za = estimate_anchor(zbr, ZERO_TARGET_PN, avoid_above_pn=step_f,
                                 min_points=RAW_ZERO_MIN_POINTS, max_reach_pn=RAW_ZERO_MAX_REACH_PN)
            if za is not None:
                rec.update(z_zero_raw_nm=round(za.value_nm, 2), zero_force_pn=round(za.force_pn, 2),
                           zero_note=za.shift_reason)
        # Full anchor: the SOP reads the fully open state on the down branch at
        # 38 pN. Here the magnet leaves a >= 38.6 pN plateau at 2 mm/s and the
        # 38 +- 4 pN window of the down move holds 1-3 samples, so the anchor is
        # taken where that same state is sampled for 10 s: the last second of
        # the plateau itself, at the plateau force.
        pl = st["plateau"]
        if st["force_pn"] >= F_OWN_STANDARD and pl.size:
            tail = pl[t[pl] >= t[pl[-1]] - END_WINDOW_S]
            e = ext[tail]
            if np.isfinite(e).sum() >= 10:
                q1, q3 = np.nanpercentile(e, [25, 75])
                rec.update(z_full_raw_nm=round(float(np.nanmedian(e)), 2),
                           full_force_pn=round(st["force_pn"], 2), full_iqr_nm=round(float(q3 - q1), 1),
                           t_full_done_s=round(float(t[pl[-1]] - t[0]), 2))
        rows.append(rec)
    return pool_zeros(rows)


def theil_sen(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    i, j = np.triu_indices(x.size, k=1)
    dx = x[j] - x[i]
    slope = float(np.median((y[j] - y[i])[dx != 0] / dx[dx != 0])) if (dx != 0).any() else 0.0
    return slope, float(np.median(y - slope * x))


def pool_zeros(rows: list[dict]) -> list[dict]:
    """Replace each step's own zero by the bead's robust zero line at that time.

    The 2.50 mm step (10 pN) never reaches 14 pN and has no own zero; it gets
    the line's value too. Outliers are dropped iteratively and flagged.
    """
    for r in rows:
        r["z_zero_own_nm"], r["zero_force_own_pn"] = r["z_zero_raw_nm"], r["zero_force_pn"]
        r["zero_outlier"], r["zero_source"] = False, "own" if np.isfinite(r["z_zero_raw_nm"]) else "line"
    own = [r for r in rows if np.isfinite(r["z_zero_own_nm"])]
    inl = list(own)
    for _ in range(5):
        if len(inl) < MIN_ZERO_STEPS:
            break
        x = np.array([r["t_plateau_start_s"] for r in inl]); y = np.array([r["z_zero_own_nm"] for r in inl])
        slope, icpt = theil_sen(x, y)
        res = {id(r): r["z_zero_own_nm"] - (slope * r["t_plateau_start_s"] + icpt) for r in own}
        new_inl = [r for r in own if abs(res[id(r)]) <= ZERO_OUTLIER_NM]
        if len(new_inl) == len(inl):
            break
        inl = new_inl
    if len(inl) < MIN_ZERO_STEPS:
        for r in rows:
            r.update(z_zero_raw_nm=np.nan, zero_force_pn=np.nan,
                     zero_note=f"fewer than {MIN_ZERO_STEPS} consistent step zeros")
        return rows
    zf = float(np.median([r["zero_force_own_pn"] for r in inl]))
    for r in rows:
        line = slope * r["t_plateau_start_s"] + icpt
        if np.isfinite(r["z_zero_own_nm"]):
            r["zero_residual_nm"] = round(r["z_zero_own_nm"] - line, 1)
            r["zero_outlier"] = abs(r["z_zero_own_nm"] - line) > ZERO_OUTLIER_NM
        r.update(z_zero_raw_nm=round(line, 2), zero_force_pn=round(zf, 2),
                 zero_line_slope_nm_per_s=round(slope, 3), n_zero_inliers=len(inl))
    return rows


def bead_reference(anc: pd.DataFrame, trk: dict) -> dict | None:
    """The highest-force step with a full anchor, tracked through it."""
    ok = anc[np.isfinite(anc.z_full_raw_nm) & np.isfinite(anc.z_zero_raw_nm)].copy()
    if ok.empty:
        return None
    bad = ok[(ok.full_iqr_nm > FULL_MAX_IQR_NM) | ok.zero_outlier.astype(bool)]
    ok = ok.drop(bad.index)
    if ok.empty:
        r = bad.loc[bad.plateau_force_pn.idxmax()]
        return {"kind": "unstable reference", "ref_step": int(r.step), "span_nm": np.nan,
                "own_span_nm": float(r.z_full_raw_nm - r.z_zero_raw_nm), "ratio": np.nan,
                "full_force_pn": float(r.full_force_pn),
                "why": (f"top plateau IQR {r.full_iqr_nm:.0f} nm > {FULL_MAX_IQR_NM:.0f}"
                        if r.full_iqr_nm > FULL_MAX_IQR_NM else
                        f"top step zero {r.zero_residual_nm:+.0f} nm off the bead's zero line")}
    lost_at = trk.get("t_loss_s", np.nan)
    if trk.get("status") == "lost":
        ok = ok[ok.t_full_done_s < lost_at]
    if ok.empty:
        return None
    r = ok.loc[ok.plateau_force_pn.idxmax()]
    span = float(r.z_full_raw_nm - r.z_zero_raw_nm)
    if not MIN_SPAN_NM <= span <= MAX_SPAN_NM:
        return {"own_span_nm": span, "ref_step": int(r.step), "kind": "unphysical span",
                "ratio": np.nan, "full_force_pn": float(r.full_force_pn), "span_nm": np.nan}
    ratio = span / expected_span_nm(float(r.zero_force_pn), float(r.full_force_pn))
    kind = ("full" if FULL_UNZIP_RATIO[0] <= ratio <= FULL_UNZIP_RATIO[1]
            else "partial" if ratio < FULL_UNZIP_RATIO[0] else "over_long")
    return {"own_span_nm": span, "span_nm": span, "ref_step": int(r.step),
            "ref_force_pn": float(r.plateau_force_pn), "full_force_pn": float(r.full_force_pn),
            "ratio": ratio, "kind": kind}


# ---------------------------------------------------------------------------
# per file
# ---------------------------------------------------------------------------
def process(d: int) -> dict:
    run = load(d)
    steps = segment_steps(run)
    trk = file_tracking(RAW / f"data_{d:03d}.h5")
    t0 = run["time_s"][0]
    anchors, refs, exts = [], {}, {}
    for roi in run["rois"]:
        tq = trk.get(roi, {"status": "unknown"})
        ext = drift_corrected_nm(run["z_um"][roi], run["time_s"])
        exts[roi] = ext
        rows = step_anchors(run, steps, roi, ext)
        a = pd.DataFrame(rows)
        a.insert(0, "bead", f"R{roi}")
        a.insert(0, "data_id", d)
        a["tracking"] = tq.get("status", "unknown")
        a["t_loss_s"] = (tq.get("t_loss_s", np.nan) - t0) if np.isfinite(tq.get("t_loss_s", np.nan)) else np.nan
        a["level_jump_um"] = tq.get("level_jump_um", 0.0)
        tq_rel = dict(tq, t_loss_s=a["t_loss_s"].iloc[0])
        ref = None
        if tq.get("status") in ("tracked", "lost") and abs(tq.get("level_jump_um", 0.0)) <= LEVEL_JUMP_UM:
            ref = bead_reference(a, tq_rel)
        refs[roi] = ref
        anchors.append(a)
    anc = pd.concat(anchors, ignore_index=True)

    # partial / over-long beads take the file's full-unzip median span
    full = [r["own_span_nm"] for r in refs.values() if r and r["kind"] == "full"]
    file_scale = float(np.median(full)) if len(full) >= MIN_FULL_BEADS_PER_FOV else np.nan
    return {"run": run, "steps": steps, "trk": trk, "anchors": anc, "refs": refs,
            "exts": exts, "file_full_scale": file_scale, "n_full": len(full)}


def finalise_scales(results: dict[int, dict]) -> float:
    full = [r["own_span_nm"] for res in results.values() for r in res["refs"].values()
            if r and r["kind"] == "full"]
    batch = float(np.median(full)) if full else np.nan
    for d, res in results.items():
        for roi, r in res["refs"].items():
            if not r or r["kind"] in ("full", "unphysical span", "unstable reference"):
                if r:
                    r["scale_note"] = (f"full unzip ({r['ratio']:.2f}x)" if r["kind"] == "full"
                                       else r.get("why", "") if r["kind"] == "unstable reference"
                                       else f"own span {r['own_span_nm']:.0f} nm outside "
                                            f"{MIN_SPAN_NM:.0f}-{MAX_SPAN_NM:.0f}")
                continue
            if np.isfinite(res["file_full_scale"]):
                r["span_nm"], scope = res["file_full_scale"], f"file median (n={res['n_full']})"
            else:
                r["span_nm"], scope = batch, f"batch median (n={len(full)})"
            r["scale_note"] = (f"{r['kind'].replace('_', '-')} (own {r['own_span_nm']:.0f} nm = "
                               f"{r['ratio']:.2f}x); scale = {scope} {r['span_nm']:.0f} nm")
    return batch


def stalls_and_samples(res: dict):
    run, steps, t = res["run"], res["steps"], res["run"]["time_s"]
    t0, f = t[0], res["run"]["force_pn"]
    d = run["data_id"]
    stall_rows, plateau_frames, trace_frames = [], [], []
    for roi in run["rois"]:
        ref = res["refs"][roi]
        a = res["anchors"][res["anchors"].bead == f"R{roi}"].set_index("step")
        if not ref or not np.isfinite(ref.get("span_nm", np.nan)):
            continue
        ext = res["exts"][roi]
        lost = a.tracking.iloc[0] == "lost"
        t_loss = a.t_loss_s.iloc[0] + t0 if lost else np.inf
        for st in steps:
            r = a.loc[st["step"]]
            if not np.isfinite(r.z_zero_raw_nm):
                continue
            # whole step cycle for the trace: up move, plateau, down move
            cyc = np.concatenate([st["up"], st["plateau"], st["down"]])
            _, luz_c = normalise(ext[cyc], f[cyc], r.z_zero_raw_nm, ref["span_nm"],
                                 r.zero_force_pn, ref["full_force_pn"])
            okc = np.isfinite(luz_c) & (t[cyc] < t_loss)
            if okc.any():
                trace_frames.append(pd.DataFrame({
                    "data_id": d, "bead": f"R{roi}", "step": st["step"],
                    "step_force_pn": round(st["force_pn"], 1),
                    "t_file_s": np.round(t[cyc][okc] - t0, 3), "force_pn": np.round(f[cyc][okc], 2),
                    "luz_bp": np.round(luz_c[okc], 1)}))
            pl = st["plateau"]
            tp = t[pl]
            _, luz = normalise(ext[pl], f[pl], r.z_zero_raw_nm, ref["span_nm"],
                               r.zero_force_pn, ref["full_force_pn"])
            keep = np.isfinite(luz) & (tp < t_loss)
            hold = float(tp[-1] - tp[0])
            t_rel = tp - tp[0]
            win_end = min(hold, STALL_AT_S)
            reached = keep.any() and t_rel[keep][-1] >= win_end - END_REACH_S
            in_win = keep & (t_rel >= win_end - END_WINDOW_S) & (t_rel <= win_end)
            end = luz[in_win] if reached else luz[:0]
            if end.size >= 3:
                q1, q3 = np.percentile(end, [25, 75])
                stall = float(np.median(end))
            else:
                q1 = q3 = stall = np.nan
            last_hold = keep & (t_rel >= hold - END_WINDOW_S)
            end_of_hold = float(np.median(luz[last_hold])) if last_hold.sum() >= 3 else np.nan
            last3 = keep & (t_rel >= win_end - 3.0) & (t_rel <= win_end)
            slope = (float(np.polyfit(t_rel[last3], luz[last3], 1)[0]) if last3.sum() >= 10 else np.nan)
            why = ("step zero off the bead's zero line (tracking level)" if bool(r.zero_outlier) else
                   "hold < 9 s" if hold < REQUIRED_PLATEAU_S else
                   "lost before the end of the window" if not reached else
                   "fork not stationary (IQR > 100 bp)" if not q3 - q1 <= MAX_END_IQR_BP else
                   "outside -150..1470 bp" if not STALL_GUARD_BP[0] <= stall <= STALL_GUARD_BP[1]
                   else "")
            stall_rows.append({
                "data_id": d, "bead": f"R{roi}", "step": st["step"], "pos_mm": st["pos_mm"],
                "force_pn": round(st["force_pn"], 1), "kind": ref["kind"],
                "span_nm": round(ref["span_nm"], 1), "zero_source": r.zero_source,
                "hold_s": round(hold, 2), "luz_start_bp": round(float(np.nanmedian(luz[keep & (t_rel <= 0.5)])), 1)
                if (keep & (t_rel <= 0.5)).any() else np.nan,
                "luz_stall_bp": round(stall, 1), "end_iqr_bp": round(float(q3 - q1), 1),
                "luz_end_of_hold_bp": round(end_of_hold, 1),
                "end_slope_bp_s": round(slope, 1),
                "still_moving": bool(abs(slope) > MOVING_SLOPE_BP_S) if np.isfinite(slope) else False,
                "used": why == "", "why_not_used": why})
            if keep.any():
                plateau_frames.append(pd.DataFrame({
                    "data_id": d, "bead": f"R{roi}", "step": st["step"],
                    "force_pn": round(st["force_pn"], 1), "t_plateau_s": np.round(t_rel[keep], 3),
                    "dt_s": np.round(frame_dt(tp)[keep], 5), "luz_bp": np.round(luz[keep], 1)}))
    cat = lambda fr: pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()
    return pd.DataFrame(stall_rows), cat(plateau_frames), cat(trace_frames)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def shade_sites(ax, horizontal: bool):
    for lo, hi in NUCLEOSOME_SITES_BP:
        (ax.axhspan if horizontal else ax.axvspan)(lo, hi, color=PALETTE["light_grey"],
                                                   alpha=0.35, lw=0, zorder=0)


def bead_figure(res: dict, bead: str, tr: pd.DataFrame, st: pd.DataFrame, pls: pd.DataFrame,
                colours: dict, title_extra: str = ""):
    """One bead: force and Luz vs time over the whole ladder, plateau dwells, stalls."""
    d, run = res["run"]["data_id"], res["run"]
    t = run["time_s"] - run["time_s"][0]
    ref = res["refs"][bead[1:]]
    fig = plt.figure(figsize=(11.0, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[3.2, 1.0], height_ratios=[1, 3], hspace=0.08, wspace=0.06)
    axf = fig.add_subplot(gs[0, 0])
    axl = fig.add_subplot(gs[1, 0], sharex=axf)
    axh = fig.add_subplot(gs[1, 1], sharey=axl)
    axf.plot(t, run["force_pn"], ".", ms=0.8, color=PALETTE["dark_grey"], rasterized=True)
    for step, g in tr.groupby("step"):
        axl.scatter(g.t_file_s, g.luz_bp, s=1.0, color=colours[float(g.step_force_pn.iloc[0])],
                    lw=0, rasterized=True)
    for s in st.itertuples():
        tt = res["steps"][s.step - 1]
        tx = float(t[tt["plateau"][0]] + min(s.hold_s, STALL_AT_S))
        if np.isfinite(s.luz_stall_bp):
            axl.plot(tx, s.luz_stall_bp, "o" if s.used else "x", ms=5,
                     mfc="none" if s.used else colours[float(s.force_pn)], mec="black", mew=0.7, zorder=5)
    bins = np.arange(-300, 1.15 * CONSTRUCT_BP + 10, 10)
    for step, g in pls.groupby("step"):
        h, e = np.histogram(g.luz_bp, bins=bins, weights=g.dt_s)
        axh.stairs(h / max(h.max(), 1e-9), e, orientation="horizontal",
                   color=colours[float(g.force_pn.iloc[0])], lw=0.9)
    for ax in (axl, axh):
        shade_sites(ax, horizontal=True)
        ax.axhline(0, color=PALETTE["light_grey"], ls="--", lw=0.6, zorder=0)
        ax.axhline(CONSTRUCT_BP, color=PALETTE["light_grey"], ls="-", lw=0.6, zorder=0)
    axl.set_ylim(-300, 1.15 * CONSTRUCT_BP)
    style_axes(axf, ylabel="Force (pN)")
    style_axes(axl, xlabel="Time (s)", ylabel="Luz (bp)")
    style_axes(axh, xlabel="Plateau dwell (norm.)")
    axf.tick_params(labelbottom=False)
    axh.tick_params(labelleft=False)
    trk = res["anchors"][res["anchors"].bead == bead].tracking.iloc[0]
    fig.suptitle(f"20261009 data_{d:03d} {bead}   tracking {trk} · {ref['scale_note']} · stall used at "
                 f"{int(st.used.sum())}/{len(st)} steps (o used, x not used; at {STALL_AT_S:.0f} s into "
                 f"the hold){title_extra}", fontsize=7, y=0.995)
    return fig


def bead_pages(res: dict, stalls: pd.DataFrame, plateau: pd.DataFrame, traces: pd.DataFrame,
               colours: dict, pdf_path: Path) -> int:
    n = 0
    with PdfPages(pdf_path) as pdf:
        for bead, tr in traces.groupby("bead", sort=True):
            fig = bead_figure(res, bead, tr, stalls[stalls.bead == bead],
                              plateau[plateau.bead == bead], colours)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)
            n += 1
    return n


# Over-long beads (own span > 1.25x a 980 bp construct) have no Luz meaning -
# a double tether or a tracking artefact - so, as on 20261004, they are kept in
# the tables and pages but left out of every statistic and summary figure.
# A force level with fewer than MIN_N_LEVEL tethers (46 pN: data_003 only) is
# listed, not tested.
STATS_KINDS = ("full", "partial")
MIN_N_LEVEL = 5


def never_unzips(st: pd.DataFrame) -> pd.Index:
    u = st[st.used]
    mx = u.groupby(["data_id", "bead"]).luz_stall_bp.max()
    return mx[mx <= NEVER_UNZIPS_BP].index


def analysed(st: pd.DataFrame) -> pd.DataFrame:
    u = st[st.used & st.kind.isin(STATS_KINDS)]
    u = u[~u.set_index(["data_id", "bead"]).index.isin(never_unzips(st))]
    n = u.groupby("force_pn").size()
    return u[u.force_pn.isin(n[n >= MIN_N_LEVEL].index)]


def summary_figures(st: pd.DataFrame, plateau: pd.DataFrame, colours: dict, levels: list[float]):
    u = analysed(st)
    # 1) stall ladder: every tether as grey points joined by a thin line,
    #    per-force median as black line
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    shade_sites(ax, horizontal=True)
    for (d, b), g in u.groupby(["data_id", "bead"]):
        g = g.sort_values("force_pn")
        ax.plot(g.force_pn, g.luz_stall_bp, "-", color=PALETTE["light_grey"], lw=0.4, zorder=1)
        ax.scatter(g.force_pn, g.luz_stall_bp, s=4, facecolors="none",
                   edgecolors=PALETTE["dark_grey"], lw=0.4, zorder=2)
    med = u.groupby("force_pn").luz_stall_bp.median()
    ax.plot(med.index, med.values, "-", color="black", lw=1.4, zorder=3)
    style_axes(ax, xlabel="Force (pN)", ylabel="Stall position Luz (bp)",
               ylim=(-200, 1.1 * CONSTRUCT_BP))
    fig.savefig(DEST / "fig_stall_ladder_all_tethers.pdf", bbox_inches="tight")
    fig.savefig(DEST / "fig_stall_ladder_all_tethers.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 2) stall histograms per force, stacked
    lv = [L for L in levels if (u.force_pn == L).any()]
    fig, axes = plt.subplots(len(lv), 1, figsize=(3.6, 0.55 * len(lv) + 0.8), sharex=True)
    axes = np.atleast_1d(axes)
    bins = np.arange(-200, 1.1 * CONSTRUCT_BP + 20, 20)
    for ax, L in zip(axes, lv):
        v = u[u.force_pn == L].luz_stall_bp
        shade_sites(ax, horizontal=False)
        ax.hist(v, bins=bins, color=colours[L], alpha=0.85, lw=0)
        ax.text(0.98, 0.78, f"{L:.1f} pN  n={len(v)}", transform=ax.transAxes, ha="right", fontsize=6.5)
        style_axes(ax, ylabel="")
        ax.set_yticks([])
    axes[-1].set_xlabel("Stall position Luz (bp)")
    fig.subplots_adjust(hspace=0.08)
    fig.savefig(DEST / "fig_stall_histograms_by_force.pdf", bbox_inches="tight")
    fig.savefig(DEST / "fig_stall_histograms_by_force.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 3) tether x force heatmap
    piv = u.pivot_table(index=["data_id", "bead"], columns="force_pn", values="luz_stall_bp")
    piv = piv.reindex(columns=lv)
    piv = piv.loc[piv.mean(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(3.6, max(2.5, 0.035 * len(piv) + 0.8)))
    im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0, vmax=CONSTRUCT_BP,
                   interpolation="nearest")
    ax.set_xticks(range(len(lv)))
    ax.set_xticklabels([f"{L:.0f}" for L in lv], fontsize=6)
    ax.set_yticks([])
    ax.set_xlabel("Force (pN)")
    ax.set_ylabel(f"Tether (n = {len(piv)})")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label("Stall Luz (bp)")
    fig.savefig(DEST / "fig_stall_tether_by_force_heatmap.pdf", bbox_inches="tight")
    fig.savefig(DEST / "fig_stall_tether_by_force_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 4) fraction past site-1 / site-2 entry / fully open, per force
    frac = u.groupby("force_pn").luz_stall_bp.agg(
        n="size",
        past_site1=lambda v: float((v > SITE1_ENTRY).mean()),
        past_site2=lambda v: float((v > SITE2_ENTRY).mean()),
        open_=lambda v: float((v >= FULL_OPEN_BP).mean()))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for col, c, lab in (("past_site1", PALETTE["red"], f"> {SITE1_ENTRY:.0f} bp (site-1 entry)"),
                        ("past_site2", PALETTE["blue"], f"> {SITE2_ENTRY:.0f} bp (site-2 entry)"),
                        ("open_", "black", f">= {FULL_OPEN_BP:.0f} bp (open)")):
        ax.plot(frac.index, frac[col], "o-", color=c, ms=3, lw=1.0, label=lab)
    style_axes(ax, xlabel="Force (pN)", ylabel="Fraction of tethers", ylim=(-0.04, 1.04))
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    fig.savefig(DEST / "fig_fraction_past_sites_by_force.pdf", bbox_inches="tight")
    fig.savefig(DEST / "fig_fraction_past_sites_by_force.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return frac


def statistics(st: pd.DataFrame, levels: list[float]) -> dict:
    from scipy.stats import spearmanr, wilcoxon
    u = analysed(st)
    piv = u.pivot_table(index=["data_id", "bead"], columns="force_pn", values="luz_stall_bp")
    lv = [L for L in levels if L in piv.columns]
    out = {"kinds": list(STATS_KINDS), "n_never_unzips_reported": int(len(never_unzips(st))),
           "n_tethers_any": int(len(piv)),
           "n_full": int(u[u.kind == "full"].drop_duplicates(["data_id", "bead"]).shape[0]),
           "n_partial": int(u[u.kind == "partial"].drop_duplicates(["data_id", "bead"]).shape[0]),
           "median_stall_bp": {f"{L:.1f}": round(float(piv[L].median()), 1) for L in lv},
           "n_per_force": {f"{L:.1f}": int(piv[L].notna().sum()) for L in lv}}
    steps = []
    for a, b in zip(lv[:-1], lv[1:]):
        pair = piv[[a, b]].dropna()
        if len(pair) >= 5:
            dlt = pair[b] - pair[a]
            p = float(wilcoxon(pair[b], pair[a]).pvalue) if (dlt != 0).any() else 1.0
            steps.append({"from_pn": a, "to_pn": b, "n": int(len(pair)),
                          "median_change_bp": round(float(dlt.median()), 1),
                          "frac_up": round(float((dlt > 0).mean()), 2),
                          "p": p, "p_bonferroni": min(1.0, p * (len(lv) - 1))})
    out["consecutive_steps"] = steps
    lo, hi = (18.0 if 18.0 in lv else lv[0]), lv[-1] if lv[-1] != 46.0 else lv[-2]
    pair = piv[[lo, hi]].dropna()
    if len(pair) >= 5:
        dlt = pair[hi] - pair[lo]
        out[f"{lo:.1f}_to_{hi:.1f}"] = {"n": int(len(pair)), "median_change_bp": round(float(dlt.median()), 1),
                                        "frac_up": round(float((dlt > 0).mean()), 2),
                                        "p": float(wilcoxon(pair[hi], pair[lo]).pvalue)}
    rho = []
    for _, row in piv.iterrows():
        v = row.dropna()
        if len(v) >= 4:
            rho.append(spearmanr(v.index.values, v.values).statistic)
    out["per_tether_spearman"] = {"n": len(rho), "median_rho": round(float(np.nanmedian(rho)), 2) if rho else None,
                                  "frac_positive": round(float(np.mean(np.array(rho) > 0)), 2) if rho else None}
    return out


# ---------------------------------------------------------------------------
def main() -> int:
    apply_style()
    DEST.mkdir(parents=True, exist_ok=True)
    titr, others = titration_files()
    print(f"titration files: {titr}\nnot analysed: {[o['data_id'] for o in others]}", flush=True)
    results = {}
    for d in titr:
        res = process(d)
        results[d] = res
        nref = sum(1 for r in res["refs"].values() if r)
        print(f"data_{d:03d}: {len(res['run']['rois'])} ROI, {len(res['steps'])} steps, "
              f"{nref} beads with a length scale, {res['n_full']} full", flush=True)
    batch_scale = finalise_scales(results)

    levels = sorted({round(st["force_pn"], 1) for res in results.values() for st in res["steps"]})
    colours = level_colours(levels)
    all_st, all_pl, all_anc, bead_rows = [], [], [], []
    for d, res in results.items():
        st, pl, tr = stalls_and_samples(res)
        all_st.append(st); all_pl.append(pl)
        all_anc.append(res["anchors"])
        for roi, r in res["refs"].items():
            bead_rows.append({"data_id": d, "bead": f"R{roi}",
                              "tracking": res["trk"].get(roi, {}).get("status", "unknown"),
                              **({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()}
                                 if r else {"kind": "no length scale"})})
        if not tr.empty:
            tr.to_csv(DEST / f"luz_traces_data_{d:03d}.csv.gz", index=False)
            n = bead_pages(res, st, pl, tr, colours,
                           DEST / f"WHZ_20261009_data_{d:03d}_ladder_by_bead.pdf")
            print(f"  data_{d:03d}: {n} bead pages", flush=True)
    st = pd.concat(all_st, ignore_index=True)
    pl = pd.concat([p for p in all_pl if not p.empty], ignore_index=True)
    anc = pd.concat(all_anc, ignore_index=True)
    beads = pd.DataFrame(bead_rows)
    anc.to_csv(DEST / "per_step_anchors.csv", index=False)
    beads.to_csv(DEST / "bead_references.csv", index=False)
    st.to_csv(DEST / "plateau_stall_by_step.csv", index=False)
    pl.to_csv(DEST / "plateau_samples.csv.gz", index=False)

    frac = summary_figures(st, pl, colours, levels)
    frac.to_csv(DEST / "fraction_past_sites_by_force.csv")
    stats = statistics(st, levels)
    (DEST / "stats.json").write_text(json.dumps(stats, indent=2))

    (DEST / "manifest.json").write_text(json.dumps({
        "script": "process_20261009.py", "raw": str(RAW), "construct": "YTrunk_2601_980bp",
        "normalised_to_bp": CONSTRUCT_BP, "nucleosome_sites_bp": NUCLEOSOME_SITES_BP,
        "calibration": f"{CALIBRATION} (1.7/1.0/0.7 mm = 18/30/46 pN), assumed unchanged for 20261009",
        "step_forces_pn": {f"{st_['pos_mm']:.2f} mm": round(st_["force_pn"], 1)
                           for st_ in next(iter(results.values()))["steps"]},
        "rest_force_pn": round(float(force_from_shift(np.array([REST_MM]), CALIBRATION)[0]), 2),
        "titration_files": titr, "not_analysed": others,
        "zero": "14 pN on each step's up move; 2.50 mm step borrows the nearest step's zero",
        "length_scale": "per bead from its highest step with a 38 pN full anchor on the down move, tracked",
        "full_unzip_ratio": FULL_UNZIP_RATIO, "batch_full_median_span_nm": round(batch_scale, 1),
        "stall": f"median Luz over {STALL_AT_S - END_WINDOW_S:.0f}-{STALL_AT_S:.0f} s of the plateau "
                 f"(IQR <= {MAX_END_IQR_BP:.0f} bp, hold >= {REQUIRED_PLATEAU_S:.0f} s, "
                 f"{STALL_GUARD_BP[0]:.0f}..{STALL_GUARD_BP[1]:.0f} bp)",
    }, indent=2, default=str))

    u = st[st.used]
    print(f"\nbeads with a length scale: {beads.kind.isin(['full', 'partial', 'over_long']).sum()} "
          f"({(beads.kind == 'full').sum()} full, {(beads.kind == 'partial').sum()} partial, "
          f"{(beads.kind == 'over_long').sum()} over-long); batch full median {batch_scale:.0f} nm")
    print(f"step stalls: {len(st)} computed, {len(u)} used; not used: "
          f"{st[~st.used].why_not_used.value_counts().to_dict()}")
    print(json.dumps(stats, indent=1)[:2500])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
