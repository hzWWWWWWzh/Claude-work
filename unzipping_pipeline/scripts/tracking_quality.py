#!/usr/bin/env python3
"""Per-ROI tracking quality, ported from the 20260913 bead-loss analysis.

A bead is judged from its own image channels, not from its force-extension
curve: the amplitude channel `A (a.u.)` collapses when the bead leaves focus or
detaches, the Z trace becomes a scatter cloud with a large rolling standard
deviation when tracking fails, and a sustained jump of more than 1.5 um is the
LUT locking onto a different focal plane. These are the criteria the user
accepted for 20260913 (`20260913_bead_loss_FEC/bead_loss_fec_20260913.py`,
thresholds in `TH`); they are reproduced here unchanged so that "tracked" means
the same thing across batches.

Why this module exists: on 20261004 the file with the highest force - the one
chosen as each bead's length reference - is exactly where beads get lost. An
anchor span measured on a scatter cloud or a LUT level jump can land anywhere
in 200-2500 nm and passes every geometric check, then rescales every other file
of that bead.
"""
from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d

TH = dict(a_ref_min=0.075, a_lost_ratio=0.60, z_roll_std_um=0.80,
          level_jump_um=1.5, lost_tail_fraction=0.95)


def rolling_std(x, w):
    m = uniform_filter1d(x, w, mode="nearest")
    m2 = uniform_filter1d(x * x, w, mode="nearest")
    return np.sqrt(np.maximum(m2 - m * m, 0))


def level_jump(t, z, i_end, fps, th=TH):
    """Largest sustained Z level switch: median of 1 s after minus 1 s before."""
    w = max(5, int(round(1.0 * fps)))
    if i_end < 3 * w:
        return 0.0
    zm = median_filter(z[:i_end], 5, mode="nearest")
    best = 0.0
    for i in range(w, i_end - w, max(1, w // 10)):
        d = float(np.median(zm[i:i + w]) - np.median(zm[i - w:i]))
        if abs(d) > abs(best):
            best = d
    return best


def file_tracking(path: Path, th=TH) -> dict[str, dict]:
    """Tracking verdict for every ROI of one HDF5 file."""
    with h5py.File(path, "r") as f:
        t = np.asarray(f["Traces/Time (s)"], float)
        names = sorted(f["Traces/Z (um)"].keys())
        Z = {n: np.asarray(f["Traces/Z (um)"][n], float) for n in names}
        A = {n: np.asarray(f["Traces/A (a.u.)"][n], float) for n in names}
    fps = (len(t) - 1) / (t[-1] - t[0])
    w2s = max(5, int(round(2.0 * fps)))
    out = {}
    for n in names:
        z = np.nan_to_num(Z[n], nan=np.nanmedian(Z[n]))
        am = median_filter(np.nan_to_num(A[n]), 21, mode="nearest")
        a_ref = float(np.percentile(am, 90))
        zsd = rolling_std(z, w2s)
        bad = (am < th["a_lost_ratio"] * a_ref) | (zsd > th["z_roll_std_um"])
        rec = {"a_ref": a_ref, "bad_fraction": float(bad.mean()),
               "max_zsd_um": float(zsd.max()), "t_loss_s": np.nan, "level_jump_um": 0.0}
        if a_ref < th["a_ref_min"]:
            rec["status"] = "no_bead"
            out[n] = rec
            continue
        tail = np.cumsum(bad[::-1])[::-1] / np.arange(len(bad), 0, -1)
        cand = np.flatnonzero(bad & (tail > th["lost_tail_fraction"]))
        i_end = len(t) if cand.size == 0 else int(cand[0])
        rec["level_jump_um"] = round(level_jump(t, z, i_end, fps, th), 2)
        if cand.size:
            rec["t_loss_s"] = round(float(t[i_end]), 2)
            rec["status"] = "lost_before_pull" if np.sum(~bad[:i_end]) < 0.5 * fps else "lost"
        elif abs(rec["level_jump_um"]) > th["level_jump_um"]:
            rec["status"] = "level_jump_artefact"
        else:
            rec["status"] = "tracked"
        out[n] = rec
    return out
