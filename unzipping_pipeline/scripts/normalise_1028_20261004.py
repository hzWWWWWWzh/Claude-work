#!/usr/bin/env python3
"""20261004: every bead, raw trace first, then normalised to 1028 (user's rule).

The user's rule (2026-10-08), for ALL beads:
  * normalise to 1028;
  * show the un-normalised trace first, then the normalised one;
  * 0 point    = the extension at 15 pN in the FIRST pull;
  * 1028 point = the maximum extension of the bead's highest-force COMPLETE
                 pull; if the bead breaks at some force, the maximum
                 extension of the pull before that force.
  The user then chose (2026-10-08) to take the 0 point in EVERY file from that
  file's own first pull (cycle 1), shared by both cycles of the file: the
  drift-corrected baselines of different files of one bead differ (C1 at 15 pN
  differs from the bead's first file by IQR -50..71 nm, > 200 nm in 14 %),
  about twice the C1-C2 difference inside a file (IQR -21..36 nm).
  normalised(file d) = 1028 * (extension - zero_d) / span,
  span = max extension - zero of the 1028-point file (one span per bead).
  If the 1028-point file's zero is > 150 nm off the bead's median zero (over
  its non-jump files, >= 3 zeros), that file uses the bead's median zero
  instead, for the span and its own normalisation (user, 2026-10-08: rule A
  with a robust zero; the zero of a last file whose molecule had not
  re-closed was raising zeros by a median 56 nm).
  The two 601 sites are drawn at 314-460 and 511-657 (SnapGene coordinates,
  the user's choice).

How the rule is made concrete here (v5, 2026-10-08; stated so it can be checked):
  * 0 point = the molecule closed but stretched as far as it goes (user: only
    the WLC stretch, no ssDNA). v5.1: one closed stretch per tether, measured
    from all its pulls (bead_zero; full rule in the comment block above
    CURVE_FROM_PN), added to each file's own 3.5-6.5 pN level (lower of C1/C2,
    stuck cycles left out). A file whose ramp shows re-zipping (extension falls
    > 100 nm while the force rises through 4-14 pN) started open and takes the
    5 pN level of the nearest earlier re-closed file. Replaces v4's per-cycle
    WLC-model zero and the 150 nm median-zero rule for the 1028 file; v5's
    150 nm "not re-closed" offset rule was dropped after verify7.
  * complete pull: an Fmax plateau tracked to its end (user: "撑过了 Fmax 平台
    就算完整"); v4 needed the file's LAST plateau. Files with a LUT level jump
    never count. A loss before the end of the last plateau is a break only if
    the bead never comes back in a later file. The 1028 point is taken from the
    highest-force file with a complete pull, up to and including the break file.
  * maximum extension: the largest running median over a +-0.5 s time window
    on that file's complete Fmax plateaus (user: the furthest unzipped length
    at the highest force), after a Hampel filter (+-0.5 s, 3 robust sigma).
  * extension: drift-corrected Z (SOP 3.1: line through the first and last
    second). For a bead lost before the end of a file the last second is
    after the loss, which tilted whole files before (median 1.7 um tail
    offset); there the drift line goes through the low-force dwells (< 1 pN)
    before the loss instead, or is a constant if those do not span 10 s.
  * forces use the batch calibration (1.7 / 1.0 / 0.7 mm = 18 / 30 / 46 pN);
    15 pN sits at 1.95 mm, below the lowest calibration point (extrapolated).

Outputs (WHZ_20261004/norm1028_v5/; v1 in norm1028/ is never
written because the user selected beads there by deleting pages):
  WHZ_20261004_norm1028_FOVxx.pdf   one page per bead: raw extension (each
      file's zero as a dashed segment with a filled triangle where it was
      taken; the 1028 point with a filled circle), normalised trace, force,
      and the plateau dwell of every file side by side
  norm1028_beads.csv                one row per bead: zero / max files and
      values, span, break, per-file verdicts, plateau-end values, `keep`
  norm1028_plateau_end.csv          per bead x file x cycle: plateau-end value
  WHZ_20261004_norm1028_v5_selected_N.pdf  the user's current selection
      (SELECTION_PDF) re-drawn with this rule, bookmarked by field and bead
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import NullLocator
from scipy.ndimage import median_filter, uniform_filter1d
from whz_plot import PALETTE, apply_style, raw_trace, style_axes

from process_tierA_sop import DRIFT_WINDOW_S, drift_corrected_nm
from process_20261004 import CALIBRATION, DEST, FOVS, NUCLEOSOME_SITES_BP, RAW
from reprocess_26_core_sop import HANDLE_BP, x_ds
from tether_stall_titration_20261004 import frame_dt, level_colours, plateaus
from unzip_raw_points import load_raw_run

# ZERO_RULE=c1_10pn (norm1028_v6): 0 = each file's C1 extension at 10 pN. Made
# on 2026-10-11 from an instruction that was meant for the 20261009 ladder data
# ("令10pN（第一个cycle）的最大extension为0点" = the 10 pN first step there), so
# the 20261004 default stays v5.3 (bead-level closed stretch).
import os as _os
ZERO_RULE = _os.environ.get("ZERO_RULE", "closed")
OUT = DEST / ("norm1028_v6" if ZERO_RULE == "c1_10pn" else "norm1028_v5_3")
# the user's latest curated selection: the merged PDF they trimmed to 66 pages
SELECTION_PDF = DEST / "WHZ_20261004_selected_75_tethers.pdf"
V1 = DEST / "norm1028"           # the user's trimmed PDFs live here; never written
MIN_RAW_RANGE_NM = 1200.0        # raw axis spans at least ~ a full unzip, so a
                                 # bead that barely moves does not look like one
SMALL_SPAN_NM = 500.0
PLATEAU_TOL_PN = 0.5
NORM_BP = 1028.0
ZERO_PN, ZERO_TOL_PN, ZERO_TOL_WIDE_PN, ZERO_MIN_SAMPLES = 15.0, 1.0, 2.0, 3
MAX_WINDOW_S = 1.0
# 0 point (user, 2026-10-08, v4): the molecule fully closed but stretched as far
# as it goes - only the handle's WLC stretch, no ssDNA. On each up ramp the
# closed state is c + HANDLE_BP * x_ds(F), with its offset c fixed from 5-10 pN
# (surely closed); a running baseline would absorb slow, gradual unzipping. The
# unzipping onset is the first sample above 10 pN that sits more than
# max(4 sigma, 30 nm) above that curve for 3 consecutive samples. The zero is
# the closed-state curve at the last closed sample (so a late onset call costs
# only the handle's small WLC stretch, ~10-20 nm); if the ramp never leaves the
# curve, at its highest force.
ONSET_FIT_FROM_PN = 5.0
ONSET_BASE_TO_PN = 10.0
ONSET_K = 3
ONSET_NSIG = 4.0
ONSET_FLOOR_NM = 30.0
ZERO_OFF_NM = 150.0
# 0 point, v5.1 (user, 2026-10-08: "其实主要问题是起点找不好"; v5 checked by
# verify7). The closed state keeps lengthening up to 9-26 pN in a bead-specific
# way (FOV07 R017: its 18 and 25.9 pN plateaus, ~80 nm above the WLC, are closed
# for the user), so the closed stretch is measured per TETHER and added to each
# file's own low-force level:
#   curve_d  = lower of C1/C2 in 1 pN force bins of the up ramp, plus the Fmax
#              plateau level (median, or its first 0.5 s if it creeps > CREEP_NM);
#   base_d   = curve_d at 3.5-6.5 pN;
#   stuck cycle (left out): C1/C2 differ by > 50 nm at 5 pN, the lower one is not
#              within 50 nm of the bead's closed level, and either (a) by 9-12.5
#              pN the gap has halved and the lower one moved up more than the
#              higher one came down, or (b) the lower one does not rise at low
#              force (< 30 nm from the ramp start to 5 pN) while the other is
#              > 120 nm higher and the bead normally rises > 100 nm;
#   not re-closed: a cycle whose extension FALLS by > REZIP_NM as the force
#              rises through 4-14 pN is re-zipping, so the file started open;
#              such a file takes the 5 pN level of the nearest earlier file
#              that re-closed (verify7: of 47 files v5 called not re-closed by a
#              150 nm offset, 27 had the shape of an offset baseline and only 10
#              showed re-zipping);
#   envelope = lowest (curve_d - base_d) over the other files; F* = the higher of
#              (ramp) the force where its 3-bin running median first leaves the
#              handle WLC stretch from 5 pN + RAMP_SLACK_NM for two bins in a row,
#              and (plateau) the highest Fmax whose plateau level is within the
#              WLC stretch + PLATEAU_SLACK_NM; closed stretch = envelope at F*;
#   zero_d   = base_d + closed stretch.
CURVE_FROM_PN = 3.0
BASE_BAND_PN = (3.5, 6.5)
BASE_PN = 5.0
CLUSTER_NM = 80.0
STUCK_GAP_NM = 50.0
STUCK_GUARD_NM = 50.0
CONVERGE_FRAC = 0.5
CONVERGE_BAND_PN = (9.0, 12.5)
FLAT_RISE_NM, FLAT_GAP_NM, BEAD_RISE_NM = 30.0, 120.0, 100.0
REZIP_NM = 130.0           # v5.3: 80-130 nm is a grey zone (own level, other choice dotted)
REZIP_GREY_NM = 80.0
STEP_FLAG_NM = 80.0
FLOOR_NM = 100.0            # a closed point (F*) is at most this far below the 5 pN level + WLC
OFFSET_FLAG_NM = 150.0      # a file this far above the bead's closed 5 pN level: shown both ways
REZIP_BAND_PN = (4.0, 14.0)
RAMP_SLACK_NM = 60.0
PLATEAU_SLACK_NM = 100.0
CLOSED_SLACK_NM = PLATEAU_SLACK_NM
CREEP_NM = 60.0
MIN_F_STAR_PN = 10.0         # closed only below this: the 0 point is not usable
# end point, v5: a pull counts as complete when its own Fmax plateau is tracked
# to its end (verify6: 9 tethers lost the bead only on C2 after a fully tracked,
# further-unzipped C1 plateau). The plateau is Hampel-filtered first (FOV07 R032:
# tracking bursts up to +500 nm set the old maximum); v5.1 also rejects a
# window maximum that stands > TRANSIENT_NM above BOTH the 1.5 s before and the
# 1.5 s after it (a burst, not a fork advance - those stay up).
HAMPEL_WINDOW_S = 1.0
HAMPEL_NSIG = 3.0
BIG_SPAN_NM = 1250.0         # above a full 1028 bp unzip at 46 pN: possible double tether
MIN_SPAN_STATS_NM = 150.0    # below this the 1028 scale multiplies noise by > 7
NOISY_PLATEAU_NM = 60.0
TRANSIENT_NM = 50.0
BURST_MIN_FMAX_PN = 20.0   # at 18 pN a 1 s excursion may be the fork itself
DRAWDOWN_NM = 150.0
# Sensitivity runs: ZERO_MODE=own - re-zipping files also keep their own 5 pN
# level (output <OUT>_ownbase); ZERO_MODE=bead - files > OFFSET_FLAG_NM above
# the bead's closed 5 pN level also take the bead's level, as in v5 (output
# <OUT>_beadbase). Which is right is not decidable from the data (verify7).
import os
ZERO_MODE = os.environ.get("ZERO_MODE", "default")
if ZERO_MODE in ("own", "bead"):
    OUT = OUT.parent / (OUT.name + {"own": "_ownbase", "bead": "_beadbase"}[ZERO_MODE])
LOW_FORCE_PN = 1.0          # low-force dwell, for the drift line of a lost bead
MIN_DWELL_SPAN_S = 10.0
END_WINDOW_S = 1.0
FILE_GAP_S = 8.0
NORM_LIM = (-400.0, 1250.0)
DWELL_BIN = 10.0
CODE = {"tracked": "T", "lost": "L", "level_jump_artefact": "J", "no_bead": "A",
        "lost_before_pull": "P"}


def denoise(z_um: np.ndarray) -> np.ndarray:
    """The SOP's smoothing (median 5 then mean 5), NaNs bridged."""
    z = np.asarray(z_um, float)
    ok = np.isfinite(z)
    if ok.sum() < 2:
        return np.full_like(z, np.nan)
    if not ok.all():
        z = np.interp(np.arange(z.size), np.flatnonzero(ok), z[ok])
    return uniform_filter1d(median_filter(z, size=5, mode="nearest"), size=5, mode="nearest")


def extension_nm(z_um, t, force, t_loss=None) -> tuple[np.ndarray, str]:
    """Drift-corrected extension. SOP line unless the bead is lost before the
    file's last second; then a line through the low-force dwells before the loss."""
    if t_loss is None or not np.isfinite(t_loss) or t_loss >= t[-1] - DRIFT_WINDOW_S:
        return drift_corrected_nm(z_um, t), "SOP first/last second"
    z = denoise(z_um)
    m = (force < LOW_FORCE_PN) & (t < t_loss - 0.5) & np.isfinite(z)
    if m.sum() >= 20 and t[m].max() - t[m].min() >= MIN_DWELL_SPAN_S:
        k, c = np.polyfit(t[m], z[m], 1)
        return (z - (k * t + c)) * 1000.0, "low-force dwells before the loss"
    head = z[t <= t[0] + DRIFT_WINDOW_S]
    return (z - np.mean(head)) * 1000.0, "constant (first second)"


def file_status(row) -> dict:
    jump = row.tracking == "level_jump_artefact" or bool(row.level_jump_bad)
    t_loss = row.t_loss_s if row.tracking == "lost" and np.isfinite(row.t_loss_s) else np.nan
    code = CODE.get(row.tracking, "?")
    if np.isfinite(t_loss):
        code += f"{t_loss:.0f}s"
    if jump and "J" not in code:
        code += "J"
    return {"tracking": row.tracking, "jump": jump, "t_loss": t_loss, "code": code,
            "present": row.tracking not in ("no_bead", "lost_before_pull")}


def last_plateau_end(run) -> float:
    """Time of the last sample at Fmax (within 0.5 pN): end of the last plateau."""
    on = np.flatnonzero(run.force_pn >= run.force_pn.max() - PLATEAU_TOL_PN)
    return float(run.time_s[on[-1]]) if on.size else float(run.time_s[-1])


def first_up(run, cycle=1):
    for c, br, seg in run.branches:
        if br == "up" and int(c) == cycle:
            return seg
    return None


def zero_point(run, ext, stat, cycle: int = 1) -> tuple[float, int, float, float]:
    """Median extension within 15 +- 1 (or 2) pN on that cycle's up ramp."""
    seg = first_up(run, cycle)
    if seg is None:
        return np.nan, 0, np.nan, np.nan
    t, f, e = run.time_s[seg], run.force_pn[seg], ext[seg]
    ok = np.isfinite(e)
    if np.isfinite(stat["t_loss"]):
        ok &= t < stat["t_loss"]
    for tol in (ZERO_TOL_PN, ZERO_TOL_WIDE_PN):
        w = ok & (np.abs(f - ZERO_PN) <= tol)
        if w.sum() >= ZERO_MIN_SAMPLES:
            return float(np.median(e[w])), int(w.sum()), tol, float(np.median(t[w]))
    return np.nan, int((ok & (np.abs(f - ZERO_PN) <= ZERO_TOL_WIDE_PN)).sum()), np.nan, np.nan


def closed_zero(run, ext, stat, cycle: int = 1):
    """Closed-state extension at the last force before unzipping starts.

    Returns (zero_nm, n_base, sigma_nm, t_zero, f_zero, f_onset) or NaNs."""
    nan = (np.nan, 0, np.nan, np.nan, np.nan, np.nan)
    seg = first_up(run, cycle)
    if seg is None:
        return nan
    t, f, e = run.time_s[seg], run.force_pn[seg], ext[seg]
    ok = np.isfinite(e) & (f >= ONSET_FIT_FROM_PN)
    if np.isfinite(stat["t_loss"]):
        ok &= t < stat["t_loss"]
    t, f, e = t[ok], f[ok], e[ok]
    base = f <= ONSET_BASE_TO_PN
    if base.sum() < 5:
        return nan
    r = e - HANDLE_BP * x_ds(f)                 # extension minus the handle's WLC stretch
    c = float(np.median(r[base]))
    sig = 1.4826 * float(np.median(np.abs(r[base] - c)))
    thr = max(ONSET_NSIG * sig, ONSET_FLOOR_NM)
    last, onset = len(f) - 1, np.nan
    above = (r - c) > thr
    for i in range(int(base.sum()), len(f)):
        if f[i] > ONSET_BASE_TO_PN and i + ONSET_K <= len(f) and above[i:i + ONSET_K].all():
            last, onset = max(i - 1, 0), float(f[i])
            break
    zero = c + HANDLE_BP * float(x_ds(np.array([f[last]]))[0])
    return zero, int(base.sum()), sig, float(t[last]), float(f[last]), onset


def fmax_plateaus(run) -> list[np.ndarray]:
    """Every contiguous run of samples within 0.5 pN of Fmax (one per cycle;
    short FOV02 holds included)."""
    idx = np.flatnonzero(run.force_pn >= run.force_pn.max() - PLATEAU_TOL_PN)
    return np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1) if idx.size else []


def complete_plateaus(run, stat) -> list[np.ndarray]:
    """Fmax plateaus the bead is still tracked at the end of (a complete pull)."""
    tl = stat["t_loss"]
    return [s for s in fmax_plateaus(run)
            if not np.isfinite(tl) or tl >= run.time_s[s[-1]]]


def hampel(t, e, window_s=HAMPEL_WINDOW_S, nsig=HAMPEL_NSIG) -> np.ndarray:
    """Samples > nsig robust sigma from their +-window/2 median are replaced by it."""
    lo = np.searchsorted(t, t - window_s / 2)
    hi = np.searchsorted(t, t + window_s / 2, side="right")
    out = e.copy()
    for j in range(e.size):
        w = e[lo[j]:hi[j]]
        m = np.median(w)
        s = 1.4826 * np.median(np.abs(w - m))
        if s > 0 and abs(e[j] - m) > nsig * s:
            out[j] = m
    return out


def plateau_noise(run, ext, segs) -> float:
    """Robust scatter (nm) about a 1 s running median on the given plateaus."""
    res = []
    for s in segs:
        t, e = run.time_s[s], ext[s]
        ok = np.isfinite(e)
        t, e = t[ok], e[ok]
        if e.size < 10:
            continue
        lo = np.searchsorted(t, t - 0.5)
        hi = np.searchsorted(t, t + 0.5, side="right")
        res.append(np.abs(e - np.array([np.median(e[a:b]) for a, b in zip(lo, hi)])))
    return float(1.4826 * np.median(np.concatenate(res))) if res else np.nan


def plateau_max(run, ext, stat, segs=None) -> tuple[float, float, float]:
    """End point (user): the furthest the fork got at the highest force - the
    largest +-0.5 s median of the extension on the file's Fmax plateau(s),
    before any loss (ramps excluded). v5: only complete plateaus (`segs`,
    default every plateau before the loss), each Hampel-filtered first."""
    t, f = run.time_s, run.force_pn
    fps = (len(t) - 1) / (t[-1] - t[0])
    if segs is None:
        ok = np.isfinite(ext) & (f >= f.max() - 0.5)
        if np.isfinite(stat["t_loss"]):
            ok &= t < stat["t_loss"]
        segs = np.split(np.flatnonzero(ok), np.flatnonzero(np.diff(np.flatnonzero(ok)) > 1) + 1)
    clean = ext.copy()
    keep = []
    for s in segs:
        s = s[np.isfinite(ext[s])]
        if np.isfinite(stat["t_loss"]):
            s = s[t[s] < stat["t_loss"]]
        if s.size >= 3:
            clean[s] = hampel(t[s], ext[s])
            keep.append(s)
    idx = np.concatenate(keep) if keep else np.array([], int)
    if not idx.size:
        return np.nan, np.nan, np.nan
    tt, ee = t[idx], clean[idx]
    lo = np.searchsorted(tt, tt - MAX_WINDOW_S / 2)
    hi = np.searchsorted(tt, tt + MAX_WINDOW_S / 2, side="right")
    need = max(3, int(0.4 * MAX_WINDOW_S * fps))
    med = np.array([np.median(ee[lo[j]:hi[j]]) if hi[j] - lo[j] >= need else -np.inf
                    for j in range(idx.size)])
    # v5.1: a maximum that stands > TRANSIENT_NM above both the 1.5 s before and
    # the 1.5 s after it is a tracking burst, not a fork advance
    masked = np.zeros(idx.size, bool)
    for j in np.argsort(-med):
        if not np.isfinite(med[j]):
            break
        if masked[j]:
            continue
        pre = ee[(tt >= tt[j] - 2.0) & (tt < tt[j] - 0.5)]
        post = ee[(tt > tt[j] + 0.5) & (tt <= tt[j] + 2.0)]
        if (f.max() > BURST_MIN_FMAX_PN and pre.size >= need and post.size >= need
                and np.median(pre) < med[j] - TRANSIENT_NM and np.median(post) < med[j] - TRANSIENT_NM):
            # v5.2: the whole burst goes, not only its top window (FOV07 R032)
            masked |= np.abs(tt - tt[j]) <= 1.0
            continue
        return float(med[j]), float(tt[j]), float(f[idx[j]])
    return np.nan, np.nan, np.nan


def max_point(run, ext, stat) -> tuple[float, float, float]:
    """Largest median over +-0.5 s of the stretched extension (F >= 15 pN,
    before any loss); a window needs at least half a second of samples."""
    t, f = run.time_s, run.force_pn
    fps = (len(t) - 1) / (t[-1] - t[0])
    ok = np.isfinite(ext) & (f >= ZERO_PN)
    if np.isfinite(stat["t_loss"]):
        ok &= t < stat["t_loss"]
    idx = np.flatnonzero(ok)
    if not idx.size:
        return np.nan, np.nan, np.nan
    tt, ee = t[idx], ext[idx]
    lo = np.searchsorted(tt, tt - MAX_WINDOW_S / 2)
    hi = np.searchsorted(tt, tt + MAX_WINDOW_S / 2, side="right")
    need = max(5, int(0.5 * MAX_WINDOW_S * fps))
    best, bi = -np.inf, -1
    for j in range(idx.size):
        if hi[j] - lo[j] >= need:
            m = float(np.median(ee[lo[j]:hi[j]]))
            if m > best:
                best, bi = m, j
    if bi < 0:
        return np.nan, np.nan, np.nan
    return best, float(tt[bi]), float(f[idx[bi]])


def cycle_curve(run, ext, stat, cycle: int) -> tuple[dict[int, float], float, float]:
    """One cycle: ({1 pN force bin: up-ramp median}, Fmax plateau level), before
    any loss. The plateau level is its median, or its first 0.5 s if it creeps
    up by > CREEP_NM (the closed state at arrival, before unzipping went on)."""
    up = first_up(run, cycle)
    if up is None:
        return {}, np.nan, np.nan
    t, f = run.time_s, run.force_pn
    fmax = float(f.max())
    good = np.isfinite(ext)
    if np.isfinite(stat["t_loss"]):
        good &= t < stat["t_loss"]
    bins = {}
    r = up[good[up] & (f[up] >= CURVE_FROM_PN) & (f[up] < fmax - PLATEAU_TOL_PN)]
    if r.size:
        b = np.round(f[r]).astype(int)
        for k in np.unique(b):
            v = ext[r][b == k]
            if v.size >= 2:
                bins[int(k)] = float(np.median(v))
    level, creep = np.nan, np.nan
    on = [s for s in fmax_plateaus(run) if s[0] >= up[0]]
    if on:
        p = on[0][good[on[0]]]
        fps = (len(t) - 1) / (t[-1] - t[0])
        hold = t[p[-1]] - t[p[0]] if p.size else 0.0
        n = max(5, int(round(min(1.0, hold / 3) * fps)))   # 1 s (v5.2), a third of short holds (v5.3)
        if p.size >= 2 * n and hold >= 0.8:
            ps, pe = np.median(ext[p][:n]), np.median(ext[p][-n:])
            creep = float(pe - ps)
            level = float(ps if creep > CREEP_NM else np.median(ext[p]))
    return bins, level, creep


def closed_band(f, slack=PLATEAU_SLACK_NM) -> np.ndarray:
    """Upper edge of the closed state above the 5 pN level: the handle's WLC
    stretch from 5 pN plus a slack."""
    f = np.asarray(f, float)
    return HANDLE_BP * (x_ds(f) - x_ds(np.array([BASE_PN]))[0]) + slack


def _band_level(bins: dict, lo: float, hi: float) -> float:
    v = [x for k, x in bins.items() if lo <= k <= hi]
    return float(np.median(v)) if v else np.nan


def _drawdown(bins: dict, lo: float, hi: float) -> float:
    """Largest fall of the extension as the force rises through lo-hi pN."""
    v = [bins[k] for k in sorted(bins) if lo <= k <= hi]
    if len(v) < 3:
        return 0.0
    v = np.array(v)
    return float(np.max(np.maximum.accumulate(v) - v))


def bead_zero(files, runs, exts, stats) -> dict | None:
    """The tether's closed, fully stretched state and every file's zero (v5.1;
    rule in the comment block above CURVE_FROM_PN)."""
    info = {}
    for d in files:
        st = stats.get(d)
        if not st or not st["present"] or st["jump"] or d not in exts:
            continue
        run, e = runs[d], exts[d]
        cy = {}
        for c in (1, 2):
            bins, level, creep = cycle_curve(run, e, st, c)
            if not bins:
                continue
            up = first_up(run, c)
            ok = np.isfinite(e[up]) & (run.force_pn[up] < 1.0)
            if np.isfinite(st["t_loss"]):
                ok &= run.time_s[up] < st["t_loss"]
            seg_v = [bins[k] for k in sorted(bins) if REZIP_BAND_PN[0] <= k <= REZIP_BAND_PN[1]]
            # noise: scatter about a 3-bin running median (a steady fall is not noise)
            if len(seg_v) >= 3:
                sv = np.array(seg_v)
                rm = np.array([np.median(sv[max(0, i - 1):i + 2]) for i in range(sv.size)])
                noise = 1.4826 * float(np.median(np.abs(sv - rm)))
            else:
                noise = 0.0
            cy[c] = {"bins": bins, "level": level, "creep": creep, "noise": noise,
                     "e0": float(np.median(e[up][ok])) if ok.sum() >= 3 else np.nan,
                     "low": _band_level(bins, *BASE_BAND_PN),
                     "mid": _band_level(bins, *CONVERGE_BAND_PN),
                     "drop": _drawdown(bins, *REZIP_BAND_PN)}
        if cy:
            info[d] = cy
    if not info:
        return None
    # provisional closed 5 pN level: the lowest file level that has another
    # within CLUSTER_NM (a lone lowest one is a shifted baseline)
    first_lv = np.sort([min(v["low"] for v in cy.values() if np.isfinite(v["low"]))
                        for cy in info.values() if any(np.isfinite(v["low"]) for v in cy.values())])
    if not first_lv.size:
        return None
    paired = ([b for b in first_lv if ((first_lv > b) & (first_lv <= b + CLUSTER_NM)).any()]
              if first_lv.size >= 3 else [])
    start = paired[0] if paired else first_lv[0]
    ref = float(np.median(first_lv[(first_lv >= start) & (first_lv <= start + CLUSTER_NM)]))
    rises = [v["low"] - v["e0"] for cy in info.values() for v in cy.values()
             if np.isfinite(v["low"]) and np.isfinite(v["e0"])]
    bead_rise = float(np.median(rises)) if rises else np.nan
    stuck, gaps, rezip, drops, flat = {}, {}, {}, {}, {}
    wlc_ = lambda f: HANDLE_BP * float(x_ds(np.array([float(f)]))[0])
    curves, plats, base, base_cyc = {}, {}, {}, {}
    for d, cy in info.items():
        if len(cy) == 2 and all(np.isfinite(cy[c]["low"]) for c in cy):
            lo_c, hi_c = sorted(cy, key=lambda c: cy[c]["low"])
            lo, hi = cy[lo_c], cy[hi_c]
            gap = hi["low"] - lo["low"]
            gaps[d] = gap
            if gap > STUCK_GAP_NM and abs(lo["low"] - ref) > STUCK_GUARD_NM:
                converging = (np.isfinite(lo["mid"]) and np.isfinite(hi["mid"])
                              and abs(hi["mid"] - lo["mid"]) < CONVERGE_FRAC * gap)
                # v5.2: the higher cycle must itself look closed (no more than
                # the WLC + ramp slack from 5 to 11 pN; FOV11 R035 D059)
                hi_closed = hi["mid"] - hi["low"] <= wlc_(11.0) - wlc_(BASE_PN) + RAMP_SLACK_NM
                if converging and hi_closed and lo["mid"] - lo["low"] > hi["low"] - hi["mid"]:
                    stuck[d] = (lo_c, "moved up to the other cycle")
                elif (np.isfinite(lo["e0"]) and lo["low"] - lo["e0"] < FLAT_RISE_NM
                      and gap > FLAT_GAP_NM and bead_rise > BEAD_RISE_NM):
                    stuck[d] = (lo_c, "no low-force rise")
        use = {c: v for c, v in cy.items() if c != stuck.get(d, (None,))[0]}
        keys = set().union(*(set(v["bins"]) for v in use.values()))
        curves[d] = {k: min(v["bins"][k] for v in use.values() if k in v["bins"]) for k in keys}
        lv = [v["level"] for v in use.values() if np.isfinite(v["level"])]
        plats[d] = min(lv) if lv else np.nan
        base[d] = _band_level(curves[d], *BASE_BAND_PN)
        base_cyc[d] = min((c for c in use if np.isfinite(use[c]["low"])),
                          key=lambda c: use[c]["low"], default=next(iter(use)))
        # a fall counts only if well above the ramp's own bin-to-bin noise (v5.2)
        drops[d] = max((v["drop"] if v["drop"] > 4 * v["noise"] else 0.0) for v in use.values())
        cr = [v["creep"] for v in use.values() if np.isfinite(v["creep"])]
        flat[d] = bool(cr) and max(abs(c_) for c_ in cr) <= CREEP_NM   # v5.3: a slide down is not flat either
    # a file started open if one of its cycles re-zips on the ramp AND its base
    # (lower cycle) sits > STUCK_GAP_NM above the nearest earlier re-closed file;
    # a closed lower cycle (FOV07 R016 D038 C1) keeps the file's own base
    # v5.2: the reference is the nearest earlier re-closed file that is not
    # itself offset (> OFFSET_FLAG_NM above the closed level), else the bead's
    # closed level - so a first file can be caught too
    clean, base_from, grey = [], {}, {}
    for d in sorted(curves):
        if not np.isfinite(base[d]):
            continue
        prevs = [x for x in clean if base[x] <= ref + OFFSET_FLAG_NM]
        if not prevs:
            # v5.3: no earlier re-closed file - keep its own level (verify9:
            # inheriting the bead level broke FOV06 R079 D026/D027)
            if drops[d] > REZIP_GREY_NM:
                grey[d] = (drops[d], None)
            clean.append(d)
            continue
        ref_d = base[prevs[-1]]
        if drops[d] > REZIP_NM and base[d] > ref_d + STUCK_GAP_NM:
            rezip[d], base_from[d] = drops[d], prevs[-1]
        else:
            if drops[d] > REZIP_GREY_NM and base[d] > ref_d + STUCK_GAP_NM:
                grey[d] = (drops[d], prevs[-1])
            clean.append(d)
    if not clean:
        return None
    env, src = {}, {}
    for d in clean:
        for k, v in curves[d].items():
            if k not in env or v - base[d] < env[k]:
                env[k], src[k] = v - base[d], d
    ks = np.array(sorted(k for k in env if k >= BASE_BAND_PN[0]), dtype=int)
    if not ks.size:
        return None
    ev = np.array([env[k] for k in ks])
    sm = np.array([np.median(ev[max(0, i - 1):i + 2]) for i in range(ev.size)])
    band = closed_band(ks, RAMP_SLACK_NM)
    wlc = lambda f: HANDLE_BP * float(x_ds(np.array([float(f)]))[0])
    # plateau levels: the lowest over the re-closed files at each Fmax
    lev = {}
    for d in clean:
        fm = round(float(runs[d].force_pn.max()), 1)
        if np.isfinite(plats[d]) and (fm not in lev or plats[d] - base[d] < lev[fm][0]):
            lev[fm] = (plats[d] - base[d], d)
    # walk up in force over ramp bins and plateaus together. A point is closed
    # if it is within its slack of the 5 pN level + WLC AND of the last closed
    # plateau + WLC (ramp bins; FOV07 R049: a slowly creeping envelope) or of
    # the last closed point + WLC (plateaus; R049's 46 pN plateau sat 140 nm
    # above its closed 39.9 pN one). The ramp stops at two bins out in a row;
    # plateaus above that still count (FOV07 R017's closed 25.9 pN plateau).
    pts = [(float(k), float(v), 0, src[int(k)]) for k, v in zip(ks, sm)]
    pts += [(fm, v[0], 1, v[1]) for fm, v in lev.items()]
    floor = lambda f: wlc(f) - wlc(BASE_PN) - FLOOR_NM     # v5.2: no closed point below this
    last_pl, last = (BASE_PN, 0.0), (BASE_PN, 0.0)
    f_star, delta, how, star_file = BASE_PN, 0.0, "none", clean[0]
    outs, stopped, closed_pl = 0, False, {}
    # v5.2: order by integer force bin, ramp before plateau, so a plateau at
    # 25.9 pN and the 26 pN ramp bin meet the same way everywhere
    for f_, v, kind, d in sorted(pts, key=lambda p: (round(p[0]), p[2])):
        if kind == 0:
            if stopped:
                continue
            ok = (v <= closed_band([f_], RAMP_SLACK_NM)[0]
                  and v <= last_pl[1] + wlc(f_) - wlc(last_pl[0]) + RAMP_SLACK_NM)
            if ok:
                outs = 0
                if v >= floor(f_):          # a dip may continue the walk, not end it
                    last = (f_, v)
                    f_star, delta, how, star_file = f_, v, "ramp", d
            else:
                outs += 1
                stopped = outs >= 2
        else:
            step = PLATEAU_SLACK_NM if not stopped else RAMP_SLACK_NM   # v5.3 (FOV07 R042, R009)
            ok = (floor(f_) <= v <= closed_band([f_], PLATEAU_SLACK_NM)[0]
                  and v <= last[1] + wlc(f_) - wlc(last[0]) + step
                  and (not stopped or flat.get(d, False)))
            if ok:
                last = last_pl = (f_, v)
                closed_pl[d] = (f_, v)
                if round(f_) >= round(f_star):
                    f_star, delta, how, star_file = f_, v, "plateau", d
    offset = {d: base[d] - ref for d in clean if base[d] > ref + OFFSET_FLAG_NM}
    zeros = {}
    for d in curves:
        run = runs[d]
        own = base[d] if np.isfinite(base[d]) else ref
        inherit = base.get(base_from.get(d), ref)
        b = own if d not in rezip or ZERO_MODE == "own" else inherit
        if d in offset and ZERO_MODE == "bead":
            b = ref
        seg = first_up(run, base_cyc[d])
        tz = run.time_s[0]
        if seg is not None:
            hit = seg[run.force_pn[seg] >= f_star - 0.5]
            tz = run.time_s[hit[0]] if hit.size else run.time_s[min(seg[-1] + 1, run.time_s.size - 1)]
        zeros[d] = (b + delta, len(clean), b, float(tz), base_cyc[d], float(f_star), d not in rezip)
    return {"zeros": zeros, "ref_base_nm": ref, "base": base, "reclosed": {d: d not in rezip for d in curves},
            "rezip": rezip, "base_from": base_from, "offset": offset, "flat": flat, "grey": grey, "f_star": float(f_star), "delta_nm": delta,
            "star_file": star_file, "star_how": how, "env_f": ks, "env_nm": ev, "env_smooth_nm": sm,
            "band_nm": band, "curves": curves, "plateaus": plats, "closed_plateaus": closed_pl,
            "stuck": {d: v[0] for d, v in stuck.items()}, "stuck_why": stuck,
            "c1c2_gap": {d: g for d, g in gaps.items() if d not in stuck}, "bead_rise_nm": bead_rise}


ZERO10_PN, ZERO10_TOL, ZERO10_TOL_WIDE = 10.0, 0.5, 1.0


def zero_c1_10pn(run, ext, stat):
    """v6 0 point: the extension where the first cycle's up ramp passes 10 pN
    (median of the samples within 10 +- 0.5 pN, +- 1 pN if fewer than 3; the
    molecule is closed there and stretched as far as it goes below unzipping).
    Falls back to C2 if C1 has no such samples. Returns the zeros tuple
    (zero, n, tol, t, cycle, force, from_c1) or None."""
    for c in (1, 2):
        seg = first_up(run, c)
        if seg is None:
            continue
        t, f, e = run.time_s[seg], run.force_pn[seg], ext[seg]
        ok = np.isfinite(e)
        if np.isfinite(stat["t_loss"]):
            ok &= t < stat["t_loss"]
        for tol in (ZERO10_TOL, ZERO10_TOL_WIDE):
            w = ok & (np.abs(f - ZERO10_PN) <= tol)
            if w.sum() >= 3:
                return (float(np.median(e[w])), int(w.sum()), tol, float(np.median(t[w])), c,
                        float(np.median(f[w])), c == 1)
    return None


def bead_scale(files, runs, exts, stats):
    """v5: one closed-state zero per tether (bead_zero; every file keeps its own
    5 pN baseline if it re-closed); one span per bead from the furthest point on
    the complete Fmax plateaus of the highest-force file with a complete pull
    (before the bead breaks)."""
    rec = {"zeros": {}, "max_file": np.nan, "zmax_nm": np.nan, "t_max": np.nan,
           "f_max_pn": np.nan, "z0_ref_nm": np.nan, "break_file": np.nan,
           "f_star": np.nan, "delta_nm": np.nan, "star_file": np.nan, "ref_base_nm": np.nan,
           "not_reclosed": "", "stuck_cycles": "", "own_base_zero": {}, "grey_rezip": "", "flags": []}
    if ZERO_RULE == "c1_10pn":
        bz = None
        for d in files:
            st = stats.get(d)
            if st and st["present"] and d in exts:
                z = zero_c1_10pn(runs[d], exts[d], st)
                if z is not None:
                    rec["zeros"][d] = z
        rec.update(f_star=ZERO10_PN, delta_nm=0.0, star_how="C1 at 10 pN")
        c2 = [f"D{d:03d}" for d, z in rec["zeros"].items() if not z[6]]
        if c2:
            rec["flags"].append("no C1 samples at 10 pN, C2 used: " + ", ".join(c2))
        # a C1 already open at 10 pN sits far above the bead's other files
        zv = np.sort([z[0] for d, z in rec["zeros"].items() if not stats[d]["jump"]])
        if zv.size >= 3:
            paired = [v for v in zv if ((zv > v) & (zv <= v + CLUSTER_NM)).any()]
            lo_ = paired[0] if paired else zv[0]
            low_level = float(np.median(zv[(zv >= lo_) & (zv <= lo_ + CLUSTER_NM)]))
            high = [f"D{d:03d} +{z[0] - low_level:.0f}" for d, z in rec["zeros"].items()
                    if z[0] - low_level > OFFSET_FLAG_NM and not stats[d]["jump"]]
            if high:
                rec["flags"].append("C1 at 10 pN > 150 nm above the bead's low files (already "
                                    "open at 10 pN?): " + ", ".join(high))
    else:
        bz = bead_zero(files, runs, exts, stats)
    if bz is not None:
        rec["zeros"] = bz["zeros"]
        rec.update(f_star=bz["f_star"], delta_nm=bz["delta_nm"], star_file=bz["star_file"],
                   star_how=bz["star_how"],
                   ref_base_nm=bz["ref_base_nm"], zero_detail=bz)
        # re-zipping files: "D036 (falls 613 nm; base of D035)"
        src_ = lambda d: (f"D{bz['base_from'][d]:03d}" if bz["base_from"][d] is not None
                          else "bead closed level")
        rec["not_reclosed"] = "; ".join(f"D{d:03d} (falls {bz['rezip'][d]:.0f} nm; base of {src_(d)})"
                                        for d in bz["rezip"])
        inh = lambda d: (bz["base"][bz["base_from"][d]] if bz["base_from"][d] is not None
                         else bz["ref_base_nm"])
        rec["stuck_cycles"] = "; ".join(f"D{d:03d}:C{c}" for d, c in bz["stuck"].items())
        big = {d: g for d, g in bz["c1c2_gap"].items() if g > 150.0 and d not in bz["rezip"]}
        rec["c1c2_gap_files"] = "; ".join(f"D{d:03d}:{g:.0f}" for d, g in big.items())
        if big:
            rec["flags"].append("C1/C2 differ > 150 nm at 5 pN (lower used): "
                                + ", ".join(f"D{d:03d}" for d in big))
        # the other choice, drawn dotted for the user: a re-zipping file's own
        # level; an offset file's bead closed level (or the reverse in the
        # sensitivity runs)
        alt = {d: (inh(d) if ZERO_MODE == "own" else bz["base"][d]) for d in bz["rezip"]}
        # grey zone (falls 80-130 nm): own level used, the inherited one dotted
        for d, (dr, frm) in bz["grey"].items():
            alt[d] = bz["base"][frm] if frm is not None else bz["ref_base_nm"]
        rec["grey_rezip"] = "; ".join(
            f"D{d:03d} (falls {dr:.0f} nm; dotted = base of "
            + (f"D{frm:03d})" if frm is not None else "bead closed level)")
            for d, (dr, frm) in bz["grey"].items())
        alt.update({d: (bz["base"][d] if ZERO_MODE == "bead" else bz["ref_base_nm"])
                    for d in bz["offset"]})
        rec["own_base_zero"] = {d: v + bz["delta_nm"] for d, v in alt.items()}
        rec["offset_files"] = "; ".join(f"D{d:03d}:+{v:.0f}" for d, v in bz["offset"].items())
        # a flat plateau below F* that reads well below 0: a step between
        # closed-looking plateaus - is the higher one closed? (user's call)
        wl = lambda f: HANDLE_BP * float(x_ds(np.array([float(f)]))[0])
        steps = []
        for d in bz["curves"]:
            pl, fm = bz["plateaus"].get(d, np.nan), float(runs[d].force_pn.max())
            if (bz["reclosed"].get(d) and bz["flat"].get(d) and np.isfinite(pl)
                    and fm < bz["f_star"] - 0.5):
                read = pl - bz["base"][d] - bz["delta_nm"] - (wl(fm) - wl(bz["f_star"]))
                if read < -STEP_FLAG_NM:
                    steps.append(f"D{d:03d} ({fm:.1f} pN) {read:.0f} nm")
        if steps:
            rec["flags"].append(f"flat plateau(s) below 0 - a step before the closed level at "
                                f"{bz['f_star']:.1f} pN (closed?): " + ", ".join(steps))
    complete, lost_in_pull = [], []
    for d in files:
        st = stats.get(d)
        if not st or not st["present"] or d not in exts:
            continue
        if np.isfinite(st["t_loss"]) and st["t_loss"] < last_plateau_end(runs[d]):
            lost_in_pull.append(d)
        # complete pull: at least one Fmax plateau tracked to its end (user,
        # 2026-10-08: "撑过了 Fmax 平台就算完整")
        if complete_plateaus(runs[d], st) and not st["jump"] and d in rec["zeros"]:
            complete.append(d)
    # a loss during a pull is a break only if the bead never comes back in a
    # later file (a bead that defocuses and is tracked again did not break)
    present_after = lambda d: any(stats.get(x, {}).get("present") for x in files if x > d)
    breaks = [d for d in lost_in_pull if not present_after(d)]
    if breaks:
        rec["break_file"] = breaks[0]
        complete = [d for d in complete if d <= breaks[0]]
    rec["complete_files"] = complete
    if complete:
        d = complete[-1]
        segs = complete_plateaus(runs[d], stats[d])
        zm, tm, fm = plateau_max(runs[d], exts[d], stats[d], segs)
        rec.update(max_file=d, zmax_nm=zm, t_max=tm, f_max_pn=fm, z0_ref_nm=rec["zeros"][d][0],
                   plateau_noise_nm=plateau_noise(runs[d], exts[d], segs))
    rec["span_nm"] = rec["zmax_nm"] - rec["z0_ref_nm"]
    rec["normalised"] = bool(np.isfinite(rec["span_nm"]) and rec["span_nm"] > 0)
    # checks for the user's eye (never acted on automatically)
    further = []
    if np.isfinite(rec["span_nm"]):
        for d in complete[:-1]:
            zm, _, _ = plateau_max(runs[d], exts[d], stats[d], complete_plateaus(runs[d], stats[d]))
            gain = zm - rec["zeros"][d][0] - rec["span_nm"]
            if np.isfinite(gain) and gain > max(50.0, 0.15 * max(rec["span_nm"], 0.0)):
                further.append(f"D{d:03d}:+{gain:.0f}")
    rec["top_not_furthest"] = "; ".join(further)
    if further:
        rec["flags"].append("a lower-force pull got further: " + ", ".join(further))
    rec["zero_unreliable"] = bool(np.isfinite(rec["f_star"]) and rec["f_star"] < MIN_F_STAR_PN)
    if rec["zero_unreliable"]:
        rec["flags"].append(f"closed only up to {rec['f_star']:.0f} pN: 0 point unusable "
                            "(left out of the statistics)")
    # reach of every complete file above its own zero (a double tether shows in
    # any file, not only the 1028 one)
    big = []
    for d in complete:
        zm, _, _ = plateau_max(runs[d], exts[d], stats[d], complete_plateaus(runs[d], stats[d]))
        if np.isfinite(zm) and zm - rec["zeros"][d][0] > BIG_SPAN_NM:
            big.append(f"D{d:03d}:{zm - rec['zeros'][d][0]:.0f}")
    rec["big_reach"] = "; ".join(big)
    if big:
        rec["flags"].append(f"reach > {BIG_SPAN_NM:.0f} nm (double tether?): " + ", ".join(big))
    if np.isfinite(rec.get("plateau_noise_nm", np.nan)) and rec["plateau_noise_nm"] > NOISY_PLATEAU_NM:
        rec["flags"].append(f"noisy top plateau ({rec['plateau_noise_nm']:.0f} nm)")
    if np.isfinite(rec["max_file"]):
        d = int(rec["max_file"])
        run, e, st = runs[d], exts[d], stats[d]
        # an incomplete plateau (bead lost on it) that got further, held >= 1 s
        starts = {int(c[0]) for c in complete_plateaus(run, st)}
        inc = [s_ for s_ in fmax_plateaus(run) if int(s_[0]) not in starts
               and np.isfinite(st["t_loss"]) and run.time_s[s_[0]] < st["t_loss"]]
        for s_ in inc:
            s_ = s_[(run.time_s[s_] < st["t_loss"]) & np.isfinite(e[s_])]
            t = run.time_s[s_]
            if s_.size and t[-1] - t[0] >= 1.0:
                late = e[s_][t >= t[-1] - 1.0]
                gain = float(np.median(late)) - rec["zmax_nm"]
                if gain > 50:
                    rec["flags"].append(f"incomplete plateau before the loss got further: +{gain:.0f} nm")
    # a fall of > DRAWDOWN_NM at constant force (not physical for a fork), on
    # any complete file (v5.2; v5.1 checked only the 1028 file)
    falls = []
    for d in complete:
        run, e, st = runs[d], exts[d], stats[d]
        for s_ in complete_plateaus(run, st):
            t, v = run.time_s[s_], e[s_]
            ok = np.isfinite(v)
            if ok.sum() < 20:
                continue
            peak, _, _ = plateau_max(run, e, st, [s_])
            dd = peak - float(np.median(v[ok][t[ok] >= t[ok][-1] - 1.0]))
            if dd > DRAWDOWN_NM:
                falls.append(f"D{d:03d}:{dd:.0f}")
                break
    if falls:
        rec["flags"].append(f"falls > {DRAWDOWN_NM:.0f} nm on a plateau: " + ", ".join(falls))
    # the break file got further before the loss (held >= 1 s)
    if (np.isfinite(rec["break_file"]) and np.isfinite(rec["span_nm"])
            and int(rec["break_file"]) in rec["zeros"] and rec["break_file"] != rec["max_file"]):
        d = int(rec["break_file"])
        run, e, st = runs[d], exts[d], stats[d]
        starts = {int(c[0]) for c in complete_plateaus(run, st)}
        for s_ in fmax_plateaus(run):
            if int(s_[0]) in starts or not np.isfinite(st["t_loss"]) or run.time_s[s_[0]] >= st["t_loss"]:
                continue
            s_ = s_[(run.time_s[s_] < st["t_loss"]) & np.isfinite(e[s_])]
            t = run.time_s[s_]
            if s_.size and t[-1] - t[0] >= 1.0:
                gain = float(np.median(e[s_][t >= t[-1] - 1.0])) - rec["zeros"][d][0] - rec["span_nm"]
                if gain > max(50.0, 0.15 * max(rec["span_nm"], 0.0)):
                    reach = gain + rec["span_nm"]
                    rec["flags"].append(
                        f"before the bead was lost, D{d:03d} reached {reach:.0f} nm (> {BIG_SPAN_NM:.0f}: "
                        "double tether or breaking?)" if reach > BIG_SPAN_NM else
                        f"before the bead was lost, D{d:03d} got +{gain:.0f} nm further")
                    break
    if np.isfinite(rec["max_file"]):
        d = int(rec["max_file"])
        # the top file never unzipped while a lower one did
        top_f = float(runs[d].force_pn.max())
        if rec["normalised"] and rec["span_nm"] < MIN_SPAN_STATS_NM and (
                further or rec["f_star"] >= top_f - 0.5):
            rec["flags"].append("top file never unzipped")
    return rec


def shade_sites(ax):
    for lo, hi in NUCLEOSOME_SITES_BP:
        ax.axhspan(lo, hi, color=PALETTE["light_grey"], alpha=0.5, lw=0, zorder=0)


def draw_page(fig, fov, bead, files, runs, exts, stats, rec, cols, norm_vals):
    gs = fig.add_gridspec(3, 2, height_ratios=[1.9, 1.9, 0.7], width_ratios=[3.6, 1.4],
                          hspace=0.12, wspace=0.05)
    ax_r = fig.add_subplot(gs[0, 0])
    ax_n = fig.add_subplot(gs[1, 0])
    ax_f = fig.add_subplot(gs[2, 0])
    ax_txt = fig.add_subplot(gs[0, 1])
    ax_h = fig.add_subplot(gs[1, 1], sharey=ax_n)
    shade_sites(ax_n)
    shade_sites(ax_h)
    offset, centres, labels, raw_pts, marks, raw_xy = 0.0, [], [], [], [], []
    dwell = []
    span = rec["span_nm"] if rec["normalised"] else np.nan
    no_zero = []
    for d in files:
        run = runs[d]
        t = run.time_s - run.time_s[0]
        f = round(float(run.force_pn.max()), 1)
        raw_trace(ax_f, t + offset, run.force_pn, color=cols[f], marker_area=0.3, alpha=0.8,
                  edge_width=0.1, filled=True, rasterized=True)
        st = stats.get(d)
        code = st["code"] if st else "?"
        star = ""
        if st and st["present"] and d in exts:
            e = exts[d]
            before = np.ones(t.size, bool) if not np.isfinite(st["t_loss"]) \
                else run.time_s < st["t_loss"]
            good = before & np.isfinite(e)
            colour = PALETTE["light_grey"] if st["jump"] else cols[f]
            raw_trace(ax_r, t[good] + offset, e[good], color=colour, marker_area=0.6, alpha=0.55,
                      edge_width=0.12, filled=not st["jump"], rasterized=True)
            if (~before).any():
                raw_trace(ax_r, t[~before] + offset, e[~before], color=PALETTE["light_grey"],
                          marker_area=0.6, alpha=0.4, edge_width=0.12, filled=True,
                          rasterized=True)
            if not st["jump"]:
                raw_pts.append(e[good])
                raw_xy.append((t[good] + offset, e[good], colour))
            z = rec["zeros"].get(d)
            if z is None:
                no_zero.append(d)
                star = "*"
            else:
                # this file's zero (v5): filled = its own 5 pN baseline + the
                # bead's closed stretch; open = the molecule had not re-closed,
                # so the bead's closed baseline is used
                zc = "black" if z[6] else PALETTE["dark_grey"]
                ax_r.hlines(z[0], offset, offset + t[-1], color=zc, lw=0.7, ls="--", zorder=4)
                ax_r.plot(z[3] - run.time_s[0] + offset, z[0], marker="^", ms=4.5,
                          mfc="black" if z[6] else "white", mec=zc, ls="none", zorder=6)
                if d in rec.get("own_base_zero", {}):
                    # the alternative: this file's own 5 pN level + closed stretch
                    ax_r.hlines(rec["own_base_zero"][d], offset, offset + t[-1],
                                color=PALETTE["grey"], lw=0.7, ls=":", zorder=4)
                if not st["jump"]:
                    marks.append(z[0])
                if np.isfinite(span):
                    n = NORM_BP * (e - z[0]) / span
                    raw_trace(ax_n, t[good] + offset, n[good], color=colour, marker_area=0.6,
                              alpha=0.55, edge_width=0.12, filled=not st["jump"],
                              rasterized=True)
                    # samples beyond the fixed axis are marked at its edge, not dropped
                    for sel, y, mk in ((good & (n > NORM_LIM[1]), NORM_LIM[1], 6),
                                       (good & (n < NORM_LIM[0]), NORM_LIM[0], 7)):
                        if sel.any():
                            ax_n.plot(t[sel] + offset, np.full(sel.sum(), y), ls="none",
                                      marker=mk, ms=1.5, color=colour, alpha=0.5,
                                      clip_on=False, rasterized=True)
                    vals = []
                    for c, seg in plateaus(run).items():
                        keep = before[seg] & np.isfinite(n[seg])
                        if keep.sum() >= 5 and not st["jump"]:
                            vals.append((n[seg][keep], frame_dt(run.time_s[seg])[keep]))
                    dwell.append((d, f, vals, st["jump"]))
            if d == rec["max_file"]:
                ax_r.hlines(rec["zmax_nm"], offset, offset + t[-1], color="black", lw=0.7,
                            ls="--", zorder=4)
                ax_r.plot(rec["t_max"] - run.time_s[0] + offset, rec["zmax_nm"], marker="o",
                          ms=5, mfc="black", mec="black", ls="none", zorder=6)
                marks.append(rec["zmax_nm"])
        centres.append(offset + t[-1] / 2)
        labels.append(f"{f:.1f}{star}\nD{d:03d} {code}")
        end = offset + t[-1]
        offset = end + FILE_GAP_S
        if d != files[-1]:
            for ax in (ax_r, ax_n, ax_f):
                ax.axvline(end + FILE_GAP_S / 2, color=PALETTE["light_grey"], lw=0.6, zorder=0)
    for y in (0.0, NORM_BP):
        ax_n.axhline(y, color="black", lw=0.7, ls="--", zorder=4)
    # side panel: plateau dwell of every file (keyed by file, coloured by force)
    edges = np.arange(NORM_LIM[0], NORM_LIM[1] + DWELL_BIN, DWELL_BIN)
    hs = []
    for d, f, vals, jump in dwell:
        h = np.zeros(edges.size - 1)
        for v, w in vals:
            inr = (v >= NORM_LIM[0]) & (v <= NORM_LIM[1])
            h += np.histogram(v[inr], bins=edges, weights=w[inr])[0]
        hs.append((d, f, h, jump, bool(vals)))
    hmax = max((h.max() for _, _, h, _, _ in hs), default=0.0) or 1.0
    mids = (edges[:-1] + edges[1:]) / 2
    for k, (d, f, h, jump, has_pl) in enumerate(hs):
        ax_h.fill_betweenx(mids, k, k + 0.85 * h / hmax, step="mid", color=cols[f], alpha=0.8, lw=0)
        vs = norm_vals.get(d, [])
        off = bool(vs) and all(not NORM_LIM[0] <= v <= NORM_LIM[1] for v in vs)
        why = ("level jump" if jump else "no plateau >= 2 s" if not has_pl
               else "out of axis" if off or not h.any() else "")
        if why:
            ax_h.text(k + 0.42, 0.5, why, rotation=90, fontsize=6,
                      transform=ax_h.get_xaxis_transform(), ha="center", va="center",
                      color=PALETTE["dark_grey"])
        for v in vs:
            if NORM_LIM[0] <= v <= NORM_LIM[1]:
                ax_h.hlines(v, k, k + 0.85, color="black", lw=1.0, zorder=4)
            else:
                # a plateau end beyond the axis: a caret at the edge, not dropped
                ax_h.plot(k + 0.42, NORM_LIM[1] if v > NORM_LIM[1] else NORM_LIM[0],
                          marker=6 if v > NORM_LIM[1] else 7, ms=4, color="black",
                          clip_on=False, zorder=6)
    ax_h.set_xlim(-0.15, max(len(hs), 1))
    ax_h.set_xticks(np.arange(len(hs)) + 0.42, [f"{f:.1f}" for _, f, _, _, _ in hs], rotation=90)
    ax_h.xaxis.set_minor_locator(NullLocator())
    # summary box
    ax_txt.axis("off")
    lines = [f"{fov} {bead}", ("0 point: each file's C1 extension at 10 pN" if ZERO_RULE == "c1_10pn"
                                else "0 point: closed, fully stretched")]
    if np.isfinite(rec["f_star"]) and ZERO_RULE != "c1_10pn":
        lines += [f"  = 5 pN baseline {rec['delta_nm']:+.0f} nm",
                  f"  (closed up to {rec['f_star']:.0f} pN, D{int(rec['star_file']):03d})"]
    if rec["not_reclosed"]:
        lines.append("  started open - extension falls on the ramp (open tri.; "
                     "dotted = own level): " + rec["not_reclosed"])
    if rec.get("grey_rezip"):
        lines.append("  may have started open (own level used): " + rec["grey_rezip"])
    if rec.get("offset_files"):
        lines.append("  5 pN level > 150 nm above the bead's closed one (dotted = "
                     "bead level): " + rec["offset_files"].replace("; ", " "))
    if rec.get("stuck_cycles"):
        lines.append("  stuck cycle left out: " + rec["stuck_cycles"].replace("; ", " "))
    if np.isfinite(rec["max_file"]):
        fm = round(float(runs[int(rec["max_file"])].force_pn.max()), 1)
        lines += [f"1028 point: data_{int(rec['max_file']):03d} ({fm:.1f} pN)",
                  f"  max {rec['zmax_nm']:.0f} nm on its plateau,",
                  f"  zero {rec['z0_ref_nm']:.0f} nm"]
    else:
        lines.append("1028 point: no complete pull with a zero")
    if np.isfinite(span):
        lines.append(f"span: {span:.0f} nm")
    if np.isfinite(rec["break_file"]):
        fb = round(float(runs[int(rec["break_file"])].force_pn.max()), 1)
        lines.append(f"lost during data_{int(rec['break_file']):03d} ({fb:.1f} pN)")
    if no_zero:
        lines.append("no zero (*): " + ", ".join(f"D{d:03d}" for d in no_zero))
    import textwrap
    lines = [w for ln in lines for w in textwrap.wrap(ln, 38, subsequent_indent="    ")]
    ax_txt.text(0.04, 0.98, "\n".join(lines), transform=ax_txt.transAxes, va="top",
                ha="left", fontsize=6.5 if len(lines) <= 13 else max(4.8, 6.5 * 13 / len(lines)))
    warn = list(rec.get("flags", []))
    if np.isfinite(span) and span < SMALL_SPAN_NM:
        warn.append(f"SMALL SPAN: {span:.0f} nm (x{NORM_BP / span:.1f})")
    if warn:
        import textwrap
        fig.text(0.01, -0.02, "\n".join(textwrap.fill("! " + w, 150) for w in warn),
                 va="top", ha="left", fontsize=6.5, color=PALETTE["red"])
    if raw_pts:
        # v5.2: per file 1-99 %; the low edge from the median file (the closed
        # baseline is alike in every file), so LUT spikes of one file
        # (FOV07 R042 D034, -3000 nm) do not set the axis - they show as carets
        q = np.array([np.nanpercentile(p_, [1.0, 99.0]) for p_ in raw_pts if np.isfinite(p_).sum() > 10])
        lo, hi = (float(np.median(q[:, 0])) - 100.0, float(q[:, 1].max())) if q.size else (0.0, 1.0)
        tops = np.sort(q[:, 1]) if q.size else np.array([])
        if tops.size >= 2 and tops[-1] - tops[-2] > 1000.0:     # FOV06 R067 D027 spikes
            hi = float(tops[-2])
        lo, hi = min([lo] + marks), max([hi] + marks)
        hi = max(hi, lo + MIN_RAW_RANGE_NM)
        pad = 0.08 * (hi - lo if hi > lo else 100)
        ax_r.set_ylim(lo - pad, hi + pad)
    for ax in (ax_r, ax_n, ax_f):
        ax.set_xlim(-1, offset - FILE_GAP_S + 1)
    # raw samples beyond the axis: carets at its edge, not dropped
    ylo, yhi = ax_r.get_ylim()
    for (tt_, ee_, col_) in raw_xy:
        for sel, y, mk in ((ee_ > yhi, yhi, 6), (ee_ < ylo, ylo, 7)):
            if sel.any():
                ax_r.plot(tt_[sel], np.full(sel.sum(), y), ls="none", marker=mk, ms=1.5,
                          color=col_, alpha=0.5, clip_on=False, rasterized=True)
    ax_n.set_ylim(*NORM_LIM)
    ax_f.set_ylim(0, max(cols) * 1.1)
    if not np.isfinite(span):
        reason = ("span <= 0 (maximum below the zero)" if np.isfinite(rec["span_nm"])
                  else "no complete pull with a 15 pN zero")
        ax_n.text(0.5, 0.5, f"not normalised: {reason}",
                  transform=ax_n.transAxes, ha="center", va="center", fontsize=8,
                  color=PALETTE["dark_grey"])
    style_axes(ax_r, ylabel="Extension (nm)")
    style_axes(ax_n, ylabel=(f"Normalised (bp)\n{NORM_BP:.0f} = {span:.0f} nm" if np.isfinite(span)
                             else "Normalised (bp)"))
    if np.isfinite(span) and span < SMALL_SPAN_NM:
        ax_n.yaxis.label.set_color(PALETTE["red"])
    style_axes(ax_f, xlabel="Time, files in acquisition order (s)", ylabel="Force (pN)")
    style_axes(ax_h, xlabel="Fmax (pN)")
    ax_r.set_xticks(centres, labels)
    ax_r.xaxis.set_minor_locator(NullLocator())
    ax_r.tick_params(axis="x", labelbottom=False, labeltop=True, labelsize=7)
    ax_n.tick_params(axis="x", labelbottom=False)
    ax_h.tick_params(labelleft=False)


def user_selection() -> dict[str, set[str]]:
    """Beads in the user's latest curated selection (SELECTION_PDF, read only)."""
    import re
    import fitz
    kept: dict[str, set[str]] = {}
    for page in fitz.open(SELECTION_PDF):
        m = re.search(r"(FOV\d\d) (R\d{3})", page.get_text())
        if m:
            kept.setdefault(m.group(1), set()).add(m.group(2))
    return kept


def user_selection_v1() -> dict[str, set[str]]:
    """Beads still in the user's trimmed v1 PDFs (norm1028/), by field. Read
    only; a field counts as trimmed when its PDF has fewer pages than beads."""
    import re
    import fitz
    kept = {}
    for pdf in sorted(V1.glob("WHZ_20261004_norm1028_FOV*.pdf")):
        doc = fitz.open(pdf)
        ids = set()
        total = None
        for page in doc:
            m = re.search(r"(FOV\d\d) (R\d{3})\s+page \d+/(\d+)", page.get_text())
            if m:
                ids.add(m.group(2))
                total = int(m.group(3))
                fov = m.group(1)
        if total is not None and doc.page_count < total:
            kept[fov] = ids
    return kept


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    apply_style()
    a = pd.read_csv(DEST / "per_cycle_anchors.csv", low_memory=False)
    a["level_jump_bad"] = a.level_jump_bad.astype(str).eq("True")
    first = a.drop_duplicates(["data_id", "roi_key"])
    levels = sorted(a.fmax_pn.round(1).unique())
    cols = level_colours(levels)
    rows, ends = [], []
    kept = user_selection()
    kept_v1 = user_selection_v1()
    print(f"your selection ({SELECTION_PDF.name}): "
          f"{ {k: len(v) for k, v in kept.items()} } = {sum(map(len, kept.values()))}", flush=True)
    for fov, files in FOVS.items():
        fa = first[first.fov == fov]
        beads = sorted(fa[~fa.tracking.isin(["no_bead", "lost_before_pull"])].bead.unique())
        if not beads:
            continue
        runs = {d: load_raw_run(RAW, d, CALIBRATION) for d in files}
        pdf_path = OUT / f"WHZ_20261004_norm1028_{fov}.pdf"
        with PdfPages(pdf_path) as pdf:
            for k, bead in enumerate(beads, start=1):
                roi = bead[1:]
                stats, exts, how = {}, {}, {}
                for d in files:
                    r = fa[(fa.data_id == d) & (fa.bead == bead)]
                    if r.empty:
                        continue
                    stats[d] = file_status(r.iloc[0])
                    if stats[d]["present"] and roi in runs[d].z_um:
                        run = runs[d]
                        tl = stats[d]["t_loss"]
                        exts[d], how[d] = extension_nm(run.z_um[roi], run.time_s, run.force_pn,
                                                       tl if np.isfinite(tl) else None)
                rec = bead_scale(files, runs, exts, stats)
                # plateau-end value of every file and cycle on the normalised scale
                norm_vals = {}
                for d in files:
                    st = stats.get(d)
                    z = rec["zeros"].get(d)
                    if not rec["normalised"] or z is None or not st or not st["present"]:
                        continue
                    run = runs[d]
                    nrm = NORM_BP * (exts[d] - z[0]) / rec["span_nm"]
                    for c, seg in plateaus(run).items():
                        t = run.time_s[seg]
                        lost_before_end = bool(np.isfinite(st["t_loss"]) and st["t_loss"] < t[-1])
                        win = t >= t[-1] - END_WINDOW_S
                        v = float(np.nanmedian(nrm[seg][win])) if not lost_before_end else np.nan
                        ends.append({"fov": fov, "bead": bead, "data_id": d,
                                     "fmax_pn": round(float(run.force_pn.max()), 1), "cycle": c,
                                     "hold_s": round(float(t[-1] - t[0]), 2),
                                     "jump_file": st["jump"], "lost_before_end": lost_before_end,
                                     "plateau_end_norm": round(v, 1)})
                        if np.isfinite(v) and not st["jump"]:
                            norm_vals.setdefault(d, []).append(v)
                fig = plt.figure(figsize=(7.09, 5.2))
                draw_page(fig, fov, bead, files, runs, exts, stats, rec, cols, norm_vals)
                fig.suptitle(f"20261004 {fov} {bead}   (bead {k} of {len(beads)} in {fov})   "
                             "(top: raw extension; middle: normalised 0-1028)", y=0.995,
                             fontsize=9)
                pdf.savefig(fig, dpi=200, bbox_inches="tight")
                plt.close(fig)
                fz = lambda d: round(float(runs[d].force_pn.max()), 1)
                row = {"fov": fov, "bead": bead, "pdf": pdf_path.name, "page": k,
                       "files_present": sum(s["present"] for s in stats.values()),
                       "verdict_by_file": "; ".join(f"D{d:03d} {fz(d):.1f}:{stats[d]['code']}"
                                                    for d in files if d in stats),
                       "zero_by_file_nm": "; ".join(f"D{d:03d}:{z[0]:.0f}"
                                                    for d, z in rec["zeros"].items()),
                       "zero_cycle_by_file": "; ".join(f"D{d:03d}:C{z[4]}"
                                                       for d, z in rec["zeros"].items()),
                       "zero_force_by_file_pn": "; ".join(f"D{d:03d}:{z[5]:.1f}"
                                                          for d, z in rec["zeros"].items()),
                       "max_file": int(rec["max_file"]) if np.isfinite(rec["max_file"]) else np.nan,
                       "max_fmax_pn": fz(int(rec["max_file"])) if np.isfinite(rec["max_file"]) else np.nan,
                       "zmax_nm": round(rec["zmax_nm"], 1),
                       "force_at_max_pn": round(rec["f_max_pn"], 1),
                       "zero_of_max_file_nm": round(rec["z0_ref_nm"], 1),
                       "span_nm": round(rec["span_nm"], 1), "normalised": rec["normalised"],
                       "small_span": bool(rec["normalised"] and rec["span_nm"] < SMALL_SPAN_NM),
                       "closed_up_to_pn": rec["f_star"],
                       "closed_stretch_nm": round(rec["delta_nm"], 1),
                       "closed_stretch_file": (int(rec["star_file"]) if np.isfinite(rec["star_file"])
                                               else np.nan),
                       "bead_base_5pn_nm": round(rec["ref_base_nm"], 1),
                       "not_reclosed_files": rec["not_reclosed"],
                       "stuck_cycles": rec.get("stuck_cycles", ""),
                       "zero_unreliable": rec.get("zero_unreliable", False),
                       "big_reach": rec.get("big_reach", ""),
                       "offset_files": rec.get("offset_files", ""),
                       "grey_rezip_files": rec.get("grey_rezip", ""),
                       "c1c2_gap_files": rec.get("c1c2_gap_files", ""),
                       "alt_zero_by_file_nm": "; ".join(f"D{d:03d}:{z:.0f}" for d, z in
                                                             rec.get("own_base_zero", {}).items()),
                       "complete_files": "; ".join(f"D{d:03d}" for d in rec.get("complete_files", [])),
                       "top_not_furthest": rec.get("top_not_furthest", ""),
                       "plateau_noise_nm": round(rec.get("plateau_noise_nm", np.nan), 1),
                       "big_span": bool(rec["normalised"] and rec["span_nm"] > BIG_SPAN_NM),
                       "closed_by": rec.get("star_how", ""),
                       "span_below_stats_min": bool(rec["normalised"]
                                                    and rec["span_nm"] < MIN_SPAN_STATS_NM),
                       "flags": " | ".join(rec.get("flags", [])),
                       "selected": "y" if bead in kept.get(fov, set()) else "",
                       "kept_in_your_v1_pdf": ("y" if bead in kept_v1.get(fov, set()) else "")
                       if fov in kept_v1 else "",
                       "break_file": int(rec["break_file"]) if np.isfinite(rec["break_file"]) else np.nan,
                       "drift_line_lost_files": "; ".join(f"D{d:03d}:{how[d]}" for d in how
                                                          if how[d] != "SOP first/last second")}
                for d, vs in norm_vals.items():
                    row[f"end_D{d:03d}_{fz(d):.1f}pN"] = round(float(np.median(vs)), 0)
                row["keep"] = ""
                rows.append(row)
        print(f"{fov}: {len(beads)} beads -> {pdf_path.name}", flush=True)
    idx = pd.DataFrame(rows)
    front = [c for c in idx.columns if not c.startswith("end_D") and c != "keep"]
    idx = idx[front + sorted(c for c in idx.columns if c.startswith("end_D")) + ["keep"]]
    idx.to_csv(OUT / "norm1028_beads.csv", index=False, encoding="utf-8-sig")
    # the selected tethers, re-drawn with this rule, in one bookmarked PDF
    # (page-by-page copy keeps it small)
    import fitz
    sel_pdf, toc, n = fitz.open(), [], 0
    for fov in sorted(kept):
        src = fitz.open(OUT / f"WHZ_20261004_norm1028_{fov}.pdf")
        rows_f = idx[(idx.fov == fov) & idx.bead.isin(kept[fov])].sort_values("page")
        toc.append([1, f"{fov} ({len(rows_f)} tethers)", n + 1])
        for r in rows_f.itertuples():
            sel_pdf.insert_pdf(src, from_page=int(r.page) - 1, to_page=int(r.page) - 1)
            n += 1
            toc.append([2, f"{r.bead}  (span {r.span_nm:.0f} nm)" if r.normalised
                        else f"{r.bead}  (not normalised)", n])
    sel_pdf.set_toc(toc)
    sel_name = f"WHZ_20261004_{OUT.name}_selected_{n}.pdf"
    sel_pdf.save(OUT / sel_name, garbage=4, deflate=True, clean=True)
    print(f"selected tethers -> {sel_name} ({n} pages)", flush=True)
    pd.DataFrame(ends).to_csv(OUT / "norm1028_plateau_end.csv", index=False)
    print(f"small span (< {SMALL_SPAN_NM:.0f} nm): {int(idx.small_span.sum())}; "
          f"{int((idx.not_reclosed_files.fillna('') != '').sum())} with a re-zipping file; "
          f"{int(idx.big_span.sum())} span > {BIG_SPAN_NM:.0f}; "
          f"{int((idx.top_not_furthest.fillna('') != '').sum())} top not furthest; 1028 point force: "
          f"{idx.force_at_max_pn.describe().round(1).to_dict()}")
    print(f"\n{len(idx)} beads; normalised {int(idx.normalised.sum())}; "
          f"no complete pull with a zero {int(idx.max_file.isna().sum())}; span <= 0 "
          f"{int((idx.span_nm <= 0).sum())}; lost during a pull {int(idx.break_file.notna().sum())}")
    print(idx.span_nm.describe().round(0).to_string())
    print(f"-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
