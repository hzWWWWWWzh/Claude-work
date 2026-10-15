#!/usr/bin/env python3
"""Raw-sample extraction for display, alongside the binned table for analysis.

The pipeline keeps two representations of the same trace and they are not
interchangeable:

  BINNED  - force bins with a median and a MAD per bin. Required by the anchor
            fit (Theil-Sen needs a per-bin value), the dWLC inversion (needs an
            uncertainty to propagate) and the peak finder (needs per-bin
            weights). Produced by process_tierA_sop.bin_branch.

  RAW     - every sample, drift-corrected and put on the cycle's own 0-890
            scale. Used for every FEC and Luz-vs-force PLOT.

Plotting the binned table was a real loss of information, not a cosmetic one.
On the 2026-09-06 batch the exponential force calibration puts adjacent samples
0.27-0.56 pN apart above 15 pN - wider than any useful bin - so the SOP's
">=2 raw points per bin" rule discarded 71% of the samples inside the peak
window, and with them a 25 pN-wide intermediate plateau that is obvious in the
raw trace (20260906 D033_R012, cycle 2). Binning is a means to a fitted number;
it was never meant to be what the reader looks at.

The normalisation is identical either way: the same per-cycle zero and scale
from the binned anchor fit are applied to the raw samples, so a raw-point FEC
and its binned counterpart sit on exactly the same axes.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from process_tierA_sop import drift_corrected_nm, force_from_shift, split_branches
from reprocess_26_core_sop import (
    HANDLE_BP, MIN_DELTA_NM_PER_BP, MIN_INVERSION_FORCE, delta_nm_per_bp, x_ds,
)

FULL_BP = 890.0


@dataclass
class RawRun:
    """One HDF5 file's raw samples, ready to normalise per ROI and cycle."""

    data_id: int
    time_s: np.ndarray
    force_pn: np.ndarray
    z_um: dict[str, np.ndarray]
    branches: list[tuple[int, str, np.ndarray]]

    @property
    def rois(self) -> list[str]:
        return sorted(self.z_um)


def load_raw_run(raw_dir: Path, data_id: int,
                 calibration: str = "sop_20260724_x1.15") -> RawRun:
    """Read one file once; callers reuse it across all of its ROIs."""
    with h5py.File(Path(raw_dir) / f"data_{data_id:03d}.h5", "r") as h5:
        tr = h5["Traces"]
        time_s = np.asarray(tr["Time (s)"], dtype=float)
        shift_mm = np.asarray(tr["Shift (mm)"], dtype=float)
        z_um = {r: np.asarray(tr["Z (um)"][r], dtype=float)
                for r in sorted(tr["Z (um)"].keys())}
    return RawRun(data_id, time_s, force_from_shift(shift_mm, calibration), z_um,
                  split_branches(time_s, shift_mm))


def zero_offsets(binned: pd.DataFrame, roi_key: str,
                 anchors: pd.DataFrame) -> dict[int, float]:
    """Recover each cycle's absolute zero from the binned table.

    The anchor table stores the span but not the zero, and the binned points
    satisfy normalized = (raw - zero) * 890/span, so the zero follows. The
    median over the cycle's bins guards against a single bad row.
    """
    spans = {int(r.cycle): r.raw_anchor_span_nm
             for r in anchors[anchors.anchor_valid].itertuples()}
    out: dict[int, float] = {}
    # Anchor tables written since 2026-10-07 carry the absolute zero; use it.
    # Recovering it from the binned table only works where bins exist, and on
    # fast ramps an ROI can have none at all.
    if "z_zero_raw_nm" in anchors:
        for r in anchors[anchors.anchor_valid].itertuples():
            if np.isfinite(r.z_zero_raw_nm):
                out[int(r.cycle)] = float(r.z_zero_raw_nm)
        if out:
            return out
    g = binned[(binned.roi_key == roi_key) & binned.anchor_valid]
    for cyc, c in g.groupby("cycle"):
        span = spans.get(int(cyc), np.nan)
        if not np.isfinite(span) or span == 0:
            continue
        out[int(cyc)] = float(np.median(
            c.extension_raw_nm - c.extension_normalized_nm * span / FULL_BP))
    return out


def normalized_raw_points(run: RawRun, roi: str, binned: pd.DataFrame,
                          anchors: pd.DataFrame) -> pd.DataFrame:
    """Every raw sample of one ROI, on that cycle's 0-890 scale.

    Returns the drift-corrected extension too, so a caller can fall back to it
    for an ROI whose anchors are not valid without re-reading the file.
    """
    roi_key = f"D{run.data_id:03d}_R{roi}"
    ext = drift_corrected_nm(run.z_um[roi], run.time_s)
    zeros = zero_offsets(binned, roi_key, anchors)
    spans = {int(r.cycle): r.raw_anchor_span_nm
             for r in anchors[anchors.anchor_valid].itertuples()}
    rows = []
    for cycle, direction, seg in run.branches:
        cyc = int(cycle)
        span, zero = spans.get(cyc), zeros.get(cyc)
        norm = ((ext[seg] - zero) * FULL_BP / span
                if span and zero is not None and np.isfinite(span)
                else np.full(seg.size, np.nan))
        rows.append(pd.DataFrame({
            "roi_key": roi_key, "cycle": cyc, "branch": direction,
            "time_s": run.time_s[seg], "force_pn": run.force_pn[seg],
            "extension_raw_nm": ext[seg], "extension_normalized_nm": norm,
        }))
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    out["normalized_valid"] = np.isfinite(out.extension_normalized_nm)
    return out


def raw_points_for_file(raw_dir: Path, data_id: int, binned: pd.DataFrame,
                        anchors: pd.DataFrame) -> pd.DataFrame:
    """All ROIs of one file, for callers that want the whole table at once."""
    run = load_raw_run(raw_dir, data_id)
    anc = anchors[anchors.data_id == data_id] if "data_id" in anchors else anchors
    frames = []
    for roi in run.rois:
        roi_key = f"D{data_id:03d}_R{roi}"
        f = normalized_raw_points(run, roi, binned, anc[anc.roi_key == roi_key])
        if not f.empty:
            frames.append(f)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not out.empty:
        out.insert(0, "data_id", data_id)
    return out


def add_dwlc_luz(points: pd.DataFrame, anchors: pd.DataFrame) -> pd.DataFrame:
    """Invert the double-WLC on raw samples (SOP section 5).

    Only the uncertainty column needs a per-bin MAD, so the inversion itself
    runs on raw samples unchanged; `luz_uncertainty_bp` is simply absent here
    and any consumer that needs it must use the binned table.
    """
    out = points.copy()
    out["luz_dwlc_bp"] = np.nan
    anc = {(r.roi_key, int(r.cycle)): r
           for r in anchors[anchors.anchor_valid].itertuples()}
    for (roi_key, cyc), idx in out.groupby(["roi_key", "cycle"]).groups.items():
        r = anc.get((roi_key, int(cyc)))
        if r is None:
            continue
        sub = out.loc[idx]
        f = sub.force_pn.values
        dx = delta_nm_per_bp(f)
        ds = HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, r.zero_force_pn)))
        e_full = (HANDLE_BP * (x_ds(np.array([r.full_force_pn]))
                               - x_ds(np.array([r.zero_force_pn])))
                  + FULL_BP * delta_nm_per_bp(np.array([r.full_force_pn])))[0]
        if not np.isfinite(e_full) or e_full == 0:
            continue
        model_scale = FULL_BP / e_full
        with np.errstate(divide="ignore", invalid="ignore"):
            luz = (sub.extension_normalized_nm.values / model_scale - ds) / dx
        bad = (f < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP) | ~np.isfinite(luz)
        out.loc[idx, "luz_dwlc_bp"] = np.where(bad, np.nan, luz)
    return out
