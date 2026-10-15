#!/usr/bin/env python3
"""Reprocess the locked 26-ROI cohort per the 20260811 SOP.

Follows `20260811_unzipnucleosomedata_complete_processing_pipeline.docx`:
force x1.15 -> branch-aware per-cycle anchoring -> 0-890 nm normalisation ->
double-WLC Luz inversion -> per-ROI panels + Luz histograms.

Two deliberate departures from the stored v2 cohort, both about the zero point:

1. The SOP fixes the zero anchor at a nominal force (14 pN up, or the recorded
   exception) and estimates it over F0 +- 0.60 pN.  Auditing that window on the
   stored cohort shows 13 of 49 ROI-cycles have a real unzipping step INSIDE it,
   up to 1357 nm and 21x the local noise.  Theil-Sen then returns a value part
   way up the step, so the normalisation origin is neither the closed state nor
   the open one.  Here the anchor force is pulled back below the first resolved
   step before estimating, which is what "zero = the still-zipped state" means.

2. D017_R034 is labelled `down12_zero` in the stored cohort, but its down branch
   only covers 33.3-42.3 pN - there is no down data at 12 pN at all, and the
   stored run silently fell back to the up branch (z_zero_branch_fallback=True).
   The SOP's own section 4.1 lists it as `up, F=12 pN`, so the SOP text is right
   and the stored label is wrong.  Treated as an up-branch anchor here.

Colour convention per the user: cycle 1 = red (NPG salmon), cycle 2 = blue
(NPG dark blue).  up = filled, down = open, points never joined by lines.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "20260731_slow_force_ramp_FEC" / "final_26_ROI_core_branch_aware_analysis" / "cohort_26_normalized_points.csv.gz"
OUT = ROOT / "20260824_core26_SOP_reprocessed"

# ---- SOP section 5: double-WLC parameters -------------------------------
KBT = 4.114           # pN nm, 25 C
DS_BP_NM = 0.34       # contour per bp
DS_LP = 50.0          # persistence length, nm
DS_S = 1000.0         # stretch modulus, pN
SS_BASE_NM = 0.58     # contour per base
SS_KUHN = 1.55        # nm
SS_S = 800.0          # pN
HANDLE_BP = 1800
FULL_BP = 890.0       # Luz = 890 bp is the fully-unzipped reference coordinate
ANCHOR_HALF_WINDOW = 0.60   # pN
MIN_ANCHOR_BINS = 3
MAX_ANCHOR_REACH_PN = 1.5   # furthest a sample may lie from the anchor force
MIN_INVERSION_FORCE = 5.0
MIN_DELTA_NM_PER_BP = 0.02

# Anchor map, SOP section 4.1.  (zero_force, zero_branch, full_force, full_branch)
DEFAULT_ANCHOR = (14.0, "up", 38.0, "down")
ANCHOR_EXCEPTIONS = {
    "D017_R034": (12.0, "up", 38.0, "down"),
    "D018_R046": (12.0, "up", 33.0, "up"),
    "D025_R026": (10.0, "down", 38.0, "down"),
}

# Colours, as specified by the user: cycle 1 red, cycle 2 blue.  This is the
# SOP 5.1 convention (C1 = NPG salmon, C2 = NPG dark blue).
C_CYCLE1 = "#E64B35"   # NPG salmon / red
C_CYCLE2 = "#3C5488"   # NPG dark blue


def x_ds(force_pn: np.ndarray) -> np.ndarray:
    """dsDNA extension per bp, nm."""
    f = np.asarray(force_pn, dtype=float)
    return DS_BP_NM * (1.0 - 0.5 * np.sqrt(KBT / (f * DS_LP)) + f / DS_S)


def _langevin(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = x[small] / 3.0
    xs = x[~small]
    out[~small] = 1.0 / np.tanh(xs) - 1.0 / xs
    return out


def x_ss(force_pn: np.ndarray) -> np.ndarray:
    """ssDNA extension per base, nm."""
    f = np.asarray(force_pn, dtype=float)
    return SS_BASE_NM * (_langevin(f * SS_KUHN / KBT) + f / SS_S)


def delta_nm_per_bp(force_pn: np.ndarray) -> np.ndarray:
    """Extension gained per unzipped bp: two ss bases minus one ds bp."""
    return 2.0 * x_ss(force_pn) - x_ds(force_pn)


def theil_sen_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Median of pairwise slopes; SOP prescribes Theil-Sen for the anchors."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        return 0.0
    i, j = np.triu_indices(x.size, k=1)
    dx = x[j] - x[i]
    ok = np.abs(dx) > 1e-9
    if not ok.any():
        return 0.0
    return float(np.median((y[j] - y[i])[ok] / dx[ok]))


@dataclass
class Anchor:
    force_pn: float
    branch: str
    value_nm: float
    n_bins: int
    slope_nm_per_pn: float
    window_pn: float
    shifted_from_pn: float
    shift_reason: str


# Below this force the tether is not yet straightened and the bin-to-bin
# scatter is large, so an apparent "step" there is elasticity, not unzipping.
MIN_STEP_SEARCH_PN = 6.0
# Only steps that could actually contaminate the anchor window matter, so the
# search runs over the window plus a short approach, not the whole branch.
STEP_SEARCH_BELOW_PN = 2.0


def first_step_force(branch: pd.DataFrame, noise_floor_nm: float,
                     target_pn: float, lag: int = 1) -> float:
    """Force of the first real step near the anchor window.

    Returns +inf when nothing threatens the window, so the nominal anchor
    stands.  Restricting the search to the neighbourhood of the target matters:
    scanning the whole branch finds low-force elasticity wobbles and drags
    anchors down to 3-6 pN, where the zero is meaningless.
    """
    lo = max(MIN_STEP_SEARCH_PN, target_pn - STEP_SEARCH_BELOW_PN - ANCHOR_HALF_WINDOW)
    hi = target_pn + ANCHOR_HALF_WINDOW
    b = branch[(branch.force_bin_pn >= lo) & (branch.force_bin_pn <= hi)].sort_values("force_bin_pn")
    e = b.extension_raw_nm.values
    f = b.force_bin_pn.values
    if e.size < 4:
        return np.inf
    # `lag` compares points that far apart. Bins are independent, so lag 1;
    # raw samples have been median5+mean5 smoothed, which spreads a jump over
    # about five samples, so callers on raw samples pass lag=5 or a 300 nm jump
    # shows up as five 60 nm steps and slips under the 80 nm floor.
    lag = max(1, min(int(lag), e.size - 1))
    steps = np.abs(e[lag:] - e[:-lag])
    # Local scatter from the same neighbourhood, so the threshold adapts to the
    # bead rather than to a global constant.
    noise = np.median(steps) if steps.size else 0.0
    thresh = max(80.0, 4.0 * max(noise, noise_floor_nm))
    hit = np.flatnonzero(steps > thresh)
    return float(f[hit[0] + lag]) if hit.size else np.inf


def estimate_anchor(branch: pd.DataFrame, target_pn: float, *,
                    avoid_above_pn: float = np.inf,
                    min_points: int = MIN_ANCHOR_BINS,
                    max_reach_pn: float = MAX_ANCHOR_REACH_PN) -> Anchor | None:
    """Theil-Sen anchor at `target_pn`, pulled below `avoid_above_pn` if needed.

    The SOP estimates over target +- 0.60 pN.  When a resolved step sits inside
    that window the estimate lands mid-step, so the target is moved down until
    the whole window clears the step.
    """
    b = branch.sort_values("force_bin_pn")
    reason, original = "", target_pn
    if np.isfinite(avoid_above_pn) and target_pn + ANCHOR_HALF_WINDOW >= avoid_above_pn:
        moved = avoid_above_pn - ANCHOR_HALF_WINDOW - 0.10
        below = b[(b.force_bin_pn < avoid_above_pn - 0.10)
                  & (b.force_bin_pn >= MIN_STEP_SEARCH_PN)]
        if len(below) >= MIN_ANCHOR_BINS:
            # keep the window inside the data that exists below the step
            moved = min(moved, float(below.force_bin_pn.max()) - ANCHOR_HALF_WINDOW)
            moved = max(moved, float(below.force_bin_pn.min()) + ANCHOR_HALF_WINDOW)
            reason = f"step at {avoid_above_pn:.2f} pN inside window; moved {original:.2f}->{moved:.2f}"
            target_pn = moved

    w = b[(b.force_bin_pn - target_pn).abs() <= ANCHOR_HALF_WINDOW]
    if len(w) < min_points:
        # Widen just enough to reach `min_points`, but never beyond
        # `max_reach_pn`. Without that cap the nearest points could sit ~10 pN
        # away when a fast ramp's bins were dropped, and the "14 pN zero" was
        # silently measured at 4 or 24 pN. The Theil-Sen slope below carries
        # the estimate back to target_pn, so a wider window costs precision,
        # not bias, as long as the state does not change inside it.
        w = b.iloc[(b.force_bin_pn - target_pn).abs().argsort()[:min_points]]
        if (len(w) < min_points
                or (w.force_bin_pn - target_pn).abs().max() > max_reach_pn):
            return None
    f = w.force_bin_pn.values
    e = w.extension_raw_nm.values
    s = theil_sen_slope(f, e)
    value = float(np.median(e - s * (f - target_pn)))
    return Anchor(target_pn, str(branch.branch.iloc[0]), value, len(w), s,
                  ANCHOR_HALF_WINDOW, original, reason)


def load() -> pd.DataFrame:
    d = pd.read_csv(SOURCE, low_memory=False)
    d = d[["roi_key", "data_id", "roi", "cycle", "branch", "force_bin_pn",
           "force_median_pn", "extension_median_um", "extension_mad_um",
           "n_points"]].copy()
    d["extension_raw_nm"] = d.extension_median_um * 1000.0
    d["extension_mad_raw_nm"] = d.extension_mad_um * 1000.0
    return d.dropna(subset=["extension_raw_nm", "force_bin_pn"])


def process() -> tuple[pd.DataFrame, pd.DataFrame]:
    d = load()
    points, anchors = [], []
    for (roi, cyc), g in d.groupby(["roi_key", "cycle"], sort=True):
        zf, zb, ff, fb = ANCHOR_EXCEPTIONS.get(roi, DEFAULT_ANCHOR)
        zbr = g[g.branch == zb]
        fbr = g[g.branch == fb]
        # An exception ROI can lack the branch its anchor names (D017_R034's
        # down branch never reaches low force); fall back and record it.
        z_fallback = f_fallback = ""
        if zbr.empty or zbr.force_bin_pn.min() > zf + ANCHOR_HALF_WINDOW:
            alt = "up" if zb == "down" else "down"
            if not g[g.branch == alt].empty:
                z_fallback = f"{zb}->{alt}"
                zb, zbr = alt, g[g.branch == alt]
        if fbr.empty:
            alt = "up" if fb == "down" else "down"
            if not g[g.branch == alt].empty:
                f_fallback = f"{fb}->{alt}"
                fb, fbr = alt, g[g.branch == alt]
        if zbr.empty or fbr.empty:
            continue

        noise_floor = float(np.nanmedian(g.extension_mad_raw_nm)) or 10.0
        step_f = first_step_force(zbr, noise_floor, zf)
        za = estimate_anchor(zbr, zf, avoid_above_pn=step_f)
        fa = estimate_anchor(fbr, ff)
        if za is None or fa is None:
            continue
        span = fa.value_nm - za.value_nm
        rec = dict(roi_key=roi, cycle=int(cyc),
                   zero_force_pn=za.force_pn, zero_branch=za.branch,
                   zero_nominal_pn=za.shifted_from_pn, zero_n_bins=za.n_bins,
                   zero_slope_nm_per_pn=round(za.slope_nm_per_pn, 3),
                   zero_shift_reason=za.shift_reason, zero_branch_fallback=z_fallback,
                   first_step_pn=(None if not np.isfinite(step_f) else round(step_f, 2)),
                   full_force_pn=fa.force_pn, full_branch=fa.branch,
                   full_n_bins=fa.n_bins, full_branch_fallback=f_fallback,
                   z_zero_raw_nm=round(za.value_nm, 2), z_full_raw_nm=round(fa.value_nm, 2),
                   raw_anchor_span_nm=round(span, 2),
                   scale_890_per_nm=(round(FULL_BP / span, 6) if span > 0 else None),
                   anchor_valid=bool(span >= 200.0))
        anchors.append(rec)
        if span <= 0:
            continue

        scale = FULL_BP / span
        p = g.copy()
        p["extension_normalized_nm"] = (p.extension_raw_nm - za.value_nm) * scale
        p["extension_mad_normalized_nm"] = p.extension_mad_raw_nm * scale

        # SOP section 5: invert the double-WLC on the normalised extension.
        f = p.force_bin_pn.values
        dx = delta_nm_per_bp(f)
        ds_term = HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, za.force_pn)))
        e_full = (HANDLE_BP * (x_ds(np.array([fa.force_pn])) - x_ds(np.array([za.force_pn])))
                  + FULL_BP * delta_nm_per_bp(np.array([fa.force_pn])))[0]
        model_scale = FULL_BP / e_full if e_full != 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            luz = (p.extension_normalized_nm.values / model_scale - ds_term) / dx
            luz_unc = np.abs(p.extension_mad_normalized_nm.values / model_scale / dx)
        bad = (f < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP) | ~np.isfinite(luz)
        luz = np.where(bad, np.nan, luz)
        luz_unc = np.where(bad, np.nan, luz_unc)
        p["luz_dwlc_bp"] = luz
        p["luz_uncertainty_bp"] = luz_unc
        p["model_scale"] = model_scale
        p["anchor_valid"] = rec["anchor_valid"]
        points.append(p)

    return pd.concat(points, ignore_index=True), pd.DataFrame(anchors)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    pts, anc = process()
    pts.to_csv(OUT / "core26_normalized_points_sop.csv.gz", index=False)
    anc.to_csv(OUT / "per_cycle_anchors.csv", index=False)

    moved = anc[anc.zero_shift_reason != ""]
    print(f"ROI x cycle processed : {len(anc)}")
    print(f"anchor_valid (span>=200 nm): {int(anc.anchor_valid.sum())}/{len(anc)}")
    print(f"zero anchors moved below a step: {len(moved)}")
    if len(moved):
        print(moved[["roi_key", "cycle", "zero_nominal_pn", "zero_force_pn",
                     "first_step_pn", "raw_anchor_span_nm"]].to_string(index=False))
    fb = anc[(anc.zero_branch_fallback != "") | (anc.full_branch_fallback != "")]
    if len(fb):
        print("\nbranch fallbacks:")
        print(fb[["roi_key", "cycle", "zero_branch_fallback", "full_branch_fallback"]].to_string(index=False))

    (OUT / "analysis_manifest.json").write_text(json.dumps({
        "sop": "20260811_unzipnucleosomedata_complete_processing_pipeline.docx",
        "source_table": str(SOURCE),
        "n_roi": int(pts.roi_key.nunique()),
        "n_roi_cycles": len(anc),
        "force_scale": 1.15,
        "full_length_bp": FULL_BP,
        "anchor_half_window_pn": ANCHOR_HALF_WINDOW,
        "min_anchor_bins": MIN_ANCHOR_BINS,
        "anchor_default": DEFAULT_ANCHOR,
        "anchor_exceptions": ANCHOR_EXCEPTIONS,
        "zero_anchors_moved_below_step": int(len(moved)),
        "dwlc": {"kBT_pN_nm": KBT, "handle_bp": HANDLE_BP, "ds_nm_per_bp": DS_BP_NM,
                 "ds_Lp_nm": DS_LP, "ds_S_pN": DS_S, "ss_nm_per_base": SS_BASE_NM,
                 "ss_kuhn_nm": SS_KUHN, "ss_S_pN": SS_S},
        "colours": {"cycle1": C_CYCLE1, "cycle2": C_CYCLE2,
                    "note": "cycle 1 red, cycle 2 blue, as specified by the user"},
    }, indent=2))
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
