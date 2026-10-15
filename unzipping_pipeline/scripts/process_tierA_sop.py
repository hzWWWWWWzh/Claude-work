#!/usr/bin/env python3
"""Tier-A candidates through the full 20260811 SOP, from raw HDF5.

The 26-ROI reprocessing worked from the stored binned table because one
existed.  Tier A spans 27 files over three dates, so this starts at the raw
TraceEditor records and follows SOP section 3 end to end:

  raw Z -> denoise -> first/last-second drift line -> force from Shift x1.15
  -> branch split on dShift/dt -> 0.10 pN bins -> per-cycle anchors
  -> 0-890 normalisation -> double-WLC Luz -> panels + Luz histograms

Anchoring, the zero-point step correction and the dWLC inversion are imported
from reprocess_26_core_sop so the two cohorts cannot drift apart.
"""
from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter, uniform_filter1d

from reprocess_26_core_sop import (
    ANCHOR_EXCEPTIONS, ANCHOR_HALF_WINDOW, DEFAULT_ANCHOR, FULL_BP, HANDLE_BP,
    MIN_ANCHOR_BINS, MIN_DELTA_NM_PER_BP, MIN_INVERSION_FORCE,
    delta_nm_per_bp, estimate_anchor, first_step_force, x_ds,
)

ROOT = Path(__file__).resolve().parent
RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/ unzipping nucleosome")
TIER_A = ROOT / "20260823_recall_screen" / "SLOW_reselected_by_core26.csv"
OUT = ROOT / "20260824_tierA_SOP_processed"

# SOP section 2: 20260724 Shift-force calibration, transferred and rescaled.
CAL_SHIFT_MM = np.array([1.7, 1.5, 1.3, 1.1, 0.9, 0.8, 0.7])
CAL_FORCE_PN = np.array([18.0, 22.0, 28.0, 33.0, 40.0, 45.0, 50.0])
FORCE_SCALE = 1.15
# SOP section 3.1
SPEED_THRESHOLD_MM_S = 0.005
MAX_GAP_S = 2.0
FORCE_BIN_PN = 0.10
MIN_POINTS_PER_BIN = 2
DRIFT_WINDOW_S = 1.0


# 20261004 (user, 2026-10-07): magnet at 1.7 / 1.0 / 0.7 mm gives 18 / 30 / 46 pN,
# with no extra scale - the SOP law above gives 20.8 / 42.3 / 57.2 pN there.
# Piecewise log-linear through these three points; outside 0.7-1.7 mm the
# nearest segment's slope is extended (no calibration point exists below 18 pN,
# so forces there - the 14 pN zero anchor, the 5 pN Luz floor - are an
# extrapolation and should be confirmed).
USER_20261004_KNOTS = ((0.7, 46.0), (1.0, 30.0), (1.7, 18.0))
CALIBRATIONS = ("sop_20260724_x1.15", "user_20261004")


def force_from_shift(shift_mm: np.ndarray, calibration: str = "sop_20260724_x1.15") -> np.ndarray:
    s = np.asarray(shift_mm, dtype=float)
    if calibration == "sop_20260724_x1.15":
        slope, intercept = np.polyfit(CAL_SHIFT_MM, np.log(CAL_FORCE_PN), 1)
        return FORCE_SCALE * np.exp(intercept + slope * s)
    if calibration != "user_20261004":
        raise ValueError(f"unknown force calibration {calibration!r}; one of {CALIBRATIONS}")
    ks = np.array([k[0] for k in USER_20261004_KNOTS])
    lf = np.log([k[1] for k in USER_20261004_KNOTS])
    lo_slope = (lf[1] - lf[0]) / (ks[1] - ks[0])
    hi_slope = (lf[-1] - lf[-2]) / (ks[-1] - ks[-2])
    out = np.interp(s, ks, lf)
    out = np.where(s < ks[0], lf[0] + lo_slope * (s - ks[0]), out)
    out = np.where(s > ks[-1], lf[-1] + hi_slope * (s - ks[-1]), out)
    return np.exp(out)


def drift_corrected_nm(z_um: np.ndarray, time_s: np.ndarray) -> np.ndarray:
    """Denoise, then remove the line joining the first and last second.

    SOP 3.1 specifies exactly this drift line, which differs from the low-force
    dwell anchoring used elsewhere in the project: it needs no dwell to exist
    and so applies uniformly across files.
    """
    z = np.asarray(z_um, dtype=float)
    finite = np.isfinite(z)
    if not finite.all():
        if finite.sum() < 2:
            return np.full_like(z, np.nan)
        z = np.interp(np.arange(z.size), np.flatnonzero(finite), z[finite])
    z = uniform_filter1d(median_filter(z, size=5, mode="nearest"), size=5, mode="nearest")

    t0, t1 = time_s[0], time_s[-1]
    head = z[time_s <= t0 + DRIFT_WINDOW_S]
    tail = z[time_s >= t1 - DRIFT_WINDOW_S]
    if head.size == 0 or tail.size == 0 or t1 <= t0:
        return (z - np.median(z)) * 1000.0
    z0, z1 = float(np.mean(head)), float(np.mean(tail))
    drift = z0 + (z1 - z0) * (time_s - t0) / (t1 - t0)
    return (z - drift) * 1000.0


def split_branches(time_s: np.ndarray, shift_mm: np.ndarray) -> list[tuple[int, str, np.ndarray]]:
    """SOP 3.1: decreasing magnet shift is increasing force; keep C1/C2."""
    speed = np.gradient(shift_mm, time_s)
    out = []
    for direction, mask in (("up", speed < -SPEED_THRESHOLD_MM_S),
                            ("down", speed > SPEED_THRESHOLD_MM_S)):
        idx = np.flatnonzero(mask)
        if idx.size == 0:
            continue
        breaks = np.flatnonzero(np.diff(time_s[idx]) > MAX_GAP_S) + 1
        segs = [s for s in np.split(idx, breaks) if s.size >= 50]
        for cycle, seg in enumerate(segs[:2], start=1):   # first two = C1, C2
            out.append((cycle, direction, seg))
    return out


def bin_branch(force: np.ndarray, ext_nm: np.ndarray,
               bin_pn: float = FORCE_BIN_PN) -> pd.DataFrame:
    """Force bins of `bin_pn`, at least two raw points each (SOP 3.1).

    The width is a parameter because the SOP's 0.10 pN assumes the 253 s ramps
    of 2026-07-31.  A faster ramp samples the same force range far more
    sparsely, and a bin that holds one point is discarded by the >=2 rule.
    """
    if force.size == 0:
        return pd.DataFrame()
    lo = np.floor(force.min() / bin_pn) * bin_pn
    hi = np.ceil(force.max() / bin_pn) * bin_pn
    edges = np.arange(lo, hi + bin_pn, bin_pn)
    which = np.digitize(force, edges) - 1
    rows = []
    for b in np.unique(which):
        if b < 0 or b >= edges.size - 1:
            continue
        m = which == b
        n = int(m.sum())
        if n < MIN_POINTS_PER_BIN:
            continue
        e = ext_nm[m]
        med = float(np.median(e))
        rows.append({"force_bin_pn": float(0.5 * (edges[b] + edges[b + 1])),
                     "force_median_pn": float(np.median(force[m])),
                     "extension_raw_nm": med,
                     "extension_mad_raw_nm": float(1.4826 * np.median(np.abs(e - med))),
                     "n_points": n})
    return pd.DataFrame(rows)


def binned_table(wanted: pd.DataFrame) -> pd.DataFrame:
    """Raw HDF5 -> per-ROI, per-cycle, per-branch binned FEC."""
    frames = []
    for (date, data_id), g in wanted.groupby(["date", "data_id"]):
        path = RAW / str(date) / f"data_{int(data_id):03d}.h5"
        if not path.exists():
            print(f"  MISSING {path}")
            continue
        with h5py.File(path, "r") as h5:
            traces = h5["Traces"]
            time_s = np.asarray(traces["Time (s)"], dtype=float)
            shift_mm = np.asarray(traces["Shift (mm)"], dtype=float)
            zgrp = traces["Z (um)"]
            force_all = force_from_shift(shift_mm)
            branches = split_branches(time_s, shift_mm)
            for roi_key, roi in zip(g.roi_key, g.roi):
                if roi not in zgrp:
                    print(f"  MISSING roi {roi} in {path.name}")
                    continue
                ext = drift_corrected_nm(np.asarray(zgrp[roi], dtype=float), time_s)
                for cycle, direction, seg in branches:
                    t = bin_branch(force_all[seg], ext[seg])
                    if t.empty:
                        continue
                    t.insert(0, "branch", direction)
                    t.insert(0, "cycle", cycle)
                    t.insert(0, "roi_key", roi_key)
                    t.insert(0, "data_id", int(data_id))
                    t.insert(0, "date", int(date))
                    frames.append(t)
        print(f"  {date}/data_{int(data_id):03d}.h5: {len(g)} ROI x {len(branches)} branches")
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def normalise_and_invert(binned: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-cycle anchoring, 0-890 normalisation and dWLC inversion."""
    points, anchors = [], []
    for (date, roi, cyc), g in binned.groupby(["date", "roi_key", "cycle"], sort=True):
        zf, zb, ff, fb = ANCHOR_EXCEPTIONS.get(roi, DEFAULT_ANCHOR) if date == 20260731 \
            else DEFAULT_ANCHOR
        zbr, fbr = g[g.branch == zb], g[g.branch == fb]
        z_fallback = f_fallback = ""
        if zbr.empty or zbr.force_bin_pn.min() > zf + ANCHOR_HALF_WINDOW:
            alt = "up" if zb == "down" else "down"
            if not g[g.branch == alt].empty:
                z_fallback, zb, zbr = f"{zb}->{alt}", alt, g[g.branch == alt]
        if fbr.empty or fbr.force_bin_pn.max() < ff - ANCHOR_HALF_WINDOW:
            alt = "up" if fb == "down" else "down"
            cand = g[g.branch == alt]
            if not cand.empty and cand.force_bin_pn.max() >= ff - ANCHOR_HALF_WINDOW:
                f_fallback, fb, fbr = f"{fb}->{alt}", alt, cand
        if zbr.empty or fbr.empty:
            continue

        noise_floor = float(np.nanmedian(g.extension_mad_raw_nm)) or 10.0
        step_f = first_step_force(zbr, noise_floor, zf)
        za = estimate_anchor(zbr, zf, avoid_above_pn=step_f)
        # The full anchor must exist in the data; clamp to the branch's reach.
        ff_eff = min(ff, float(fbr.force_bin_pn.max()) - ANCHOR_HALF_WINDOW) \
            if fbr.force_bin_pn.max() < ff + ANCHOR_HALF_WINDOW else ff
        fa = estimate_anchor(fbr, ff_eff)
        if za is None or fa is None:
            continue
        span = fa.value_nm - za.value_nm
        anchors.append(dict(date=int(date), roi_key=roi, cycle=int(cyc),
            zero_force_pn=za.force_pn, zero_branch=za.branch,
            zero_nominal_pn=za.shifted_from_pn, zero_n_bins=za.n_bins,
            zero_shift_reason=za.shift_reason, zero_branch_fallback=z_fallback,
            first_step_pn=(None if not np.isfinite(step_f) else round(step_f, 2)),
            full_force_pn=fa.force_pn, full_branch=fa.branch, full_n_bins=fa.n_bins,
            full_force_nominal_pn=ff, full_branch_fallback=f_fallback,
            z_zero_raw_nm=round(za.value_nm, 2), z_full_raw_nm=round(fa.value_nm, 2),
            raw_anchor_span_nm=round(span, 2),
            scale_890_per_nm=(round(FULL_BP / span, 6) if span > 0 else None),
            anchor_valid=bool(span >= 200.0)))
        if span <= 0:
            continue

        scale = FULL_BP / span
        p = g.copy()
        p["extension_normalized_nm"] = (p.extension_raw_nm - za.value_nm) * scale
        p["extension_mad_normalized_nm"] = p.extension_mad_raw_nm * scale
        f = p.force_bin_pn.values
        dx = delta_nm_per_bp(f)
        ds_term = HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, za.force_pn)))
        e_full = (HANDLE_BP * (x_ds(np.array([fa.force_pn])) - x_ds(np.array([za.force_pn])))
                  + FULL_BP * delta_nm_per_bp(np.array([fa.force_pn])))[0]
        model_scale = FULL_BP / e_full if e_full != 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            luz = (p.extension_normalized_nm.values / model_scale - ds_term) / dx
            unc = np.abs(p.extension_mad_normalized_nm.values / model_scale / dx)
        bad = (f < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP) | ~np.isfinite(luz)
        p["luz_dwlc_bp"] = np.where(bad, np.nan, luz)
        p["luz_uncertainty_bp"] = np.where(bad, np.nan, unc)
        p["anchor_valid"] = span >= 200.0
        points.append(p)
    return (pd.concat(points, ignore_index=True) if points else pd.DataFrame(),
            pd.DataFrame(anchors))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    tier = pd.read_csv(TIER_A)
    a = tier[tier.band.str.startswith("A")].copy()
    a["data_id"] = a.roi_key.str.extract(r"D(\d+)_")[0].astype(int)
    a["roi"] = a.roi_key.str.extract(r"_R(\d+)")[0]
    print(f"Tier A ROIs: {len(a)} over {a.date.nunique()} dates, "
          f"{a.groupby(['date','data_id']).ngroups} files")

    binned = binned_table(a[["date", "data_id", "roi_key", "roi"]])
    print(f"binned rows: {len(binned)}")
    binned.to_csv(OUT / "tierA_binned_FEC.csv.gz", index=False)

    pts, anc = normalise_and_invert(binned)
    pts = pts.merge(a[["date", "roi_key", "user_core_26", "score"]],
                    on=["date", "roi_key"], how="left")
    anc = anc.merge(a[["date", "roi_key", "user_core_26", "score"]],
                    on=["date", "roi_key"], how="left")
    pts.to_csv(OUT / "tierA_normalized_points_sop.csv.gz", index=False)
    anc.to_csv(OUT / "per_cycle_anchors.csv", index=False)

    moved = anc[anc.zero_shift_reason.fillna("") != ""]
    print(f"\nROI x cycle : {len(anc)}   anchor_valid: {int(anc.anchor_valid.sum())}")
    print(f"ROIs with at least one valid cycle: {anc[anc.anchor_valid].roi_key.nunique()}")
    print(f"zero anchors moved below a step   : {len(moved)}")
    print(f"branch fallbacks (zero / full)    : "
          f"{int((anc.zero_branch_fallback.fillna('')!='').sum())} / "
          f"{int((anc.full_branch_fallback.fillna('')!='').sum())}")

    (OUT / "analysis_manifest.json").write_text(json.dumps({
        "sop": "20260811_unzipnucleosomedata_complete_processing_pipeline.docx",
        "cohort": "Tier A of 20260823_recall_screen/SLOW_reselected_by_core26.csv",
        "n_roi_requested": len(a), "n_roi_cycles": len(anc),
        "raw_source": str(RAW), "force_scale": FORCE_SCALE,
        "force_calibration_shift_mm": CAL_SHIFT_MM.tolist(),
        "force_calibration_force_pn": CAL_FORCE_PN.tolist(),
        "force_bin_pn": FORCE_BIN_PN, "min_points_per_bin": MIN_POINTS_PER_BIN,
        "branch_speed_threshold_mm_s": SPEED_THRESHOLD_MM_S,
        "drift_window_s": DRIFT_WINDOW_S,
        "full_length_bp": FULL_BP,
        "anchor_default": DEFAULT_ANCHOR, "anchor_exceptions": ANCHOR_EXCEPTIONS,
        "zero_anchors_moved_below_step": int(len(moved)),
        "colours": {"cycle1": "#E64B35 red", "cycle2": "#3C5488 blue"},
    }, indent=2))
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
