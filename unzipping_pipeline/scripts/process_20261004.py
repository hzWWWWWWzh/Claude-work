#!/usr/bin/env python3
"""20261004 force-titration batch: SOP processing, raw-point WHZ panels.

This batch is a titration, not a set of independent files. Each field of view
was pulled several times with the maximum force stepped up (18.0 -> 25.9 ->
27.9 -> 30.0 -> 34.6 -> 39.9 -> 46.0 pN on the batch calibration; 20.8 -> 57.2
pN on the old SOP law), and the ROI positions are identical
across the files of a field (0.0 px), so ROI R### in data_009 and in data_014 is
the same bead.

That changes how every file must be scaled. Normalising each file, or even
each cycle, on its own 38 pN full anchor gives one bead several different
lengths (spans of 268-428 nm for one bead in the first run), and a file whose
maximum force is low may never unzip fully - at 20 pN the fork can be stuck
before the nucleosome. Normalising such a file on its own anchor stretches a
partial opening onto the full construct length.

Rule: each bead gets ONE length scale, taken from the file in its field of view
with the HIGHEST maximum force in which the bead still has a valid own anchor
with a physical span (200-2500 nm; see MAX_SPAN_NM)
(the down branch passing 38 pN after reaching that Fmax is fully open). The
span is the median of that file's valid cycles. Every file of the bead - the
reference file included - is normalised with that span; the zero always comes
from the file and cycle itself (14 pN up, pulled below any step in its window).
Only a cycle in which the bead is TRACKED through its anchors can serve as the
reference (tracking_quality.py, the 20260913 bead-loss criteria; a bead lost
after the down branch has left the 38 pN window still counts, since on these
ramps most losses happen in the final near-zero-force dwell): the highest-force file is where
beads get lost, and a span measured on a scatter cloud or a LUT level jump
passes every geometric check. A file in which the bead is absent or was lost
before the pull is not normalised at all.
Forces use the batch calibration the user gave on 2026-10-07 (1.7 / 1.0 / 0.7 mm
= 18 / 30 / 46 pN; see CALIBRATION), not the SOP law x1.15, which overstated
them by 16-41 % (41 % at 1.0 mm, 24 % at 0.7 mm). Fields whose highest pull stays below 38.6 pN on
this scale (FOV04, FOV05: 34.6 pN) have no full anchor and are not normalised.
The construct is YTrunk_2601_980bp, normalised to 980 bp, with two 601 sites at
314-460 and 511-657 (shaded on the Luz panels). A fully unzipped 980 bp
construct spans 797 nm between the anchors (dWLC). A bead whose reference span
is below 0.8x that never unzipped fully even at its highest force - its own
span measures only how far the fork got - so it is not a 980 bp scale; it takes
the median span of the fully unzipped beads in its field (the user's choice,
2026-10-07), or of the whole batch if the field has fewer than three.
A bead with no file reaching 38.6 pN with a valid anchor is not normalised: its
pages show drift-corrected extension. Nothing is borrowed across fields of view.

The tracking verdict gates every file, not only the reference (fixed
2026-10-07: before, a non-reference file was dropped only when the bead was
absent, and 251 of 1862 normalised cycles ended up with a Luz at Fmax outside
-150..1274 bp, down to -9558 and up to 22468). A file with a LUT level jump
(> 1.5 um before any loss) is not normalised; a lost bead is normalised only up
to the loss, and a cycle whose zero window (14 pN + 2.5) comes after the loss is
dropped. Last, a cycle whose Luz at its highest force (median, Fmax - 3 pN) is
outside LUZ_GUARD_BP cannot be one tether of this construct and is dropped; the
upper bound is not applied to over-long beads, whose Luz has no meaning anyway.
`anchor_source` names the reference file for every cycle, and `own_span_nm`
keeps the file's own span for audit.

Display uses raw samples (PIPELINE_raw_vs_binned_zh.md); the 0.25 pN bins only
feed the anchor fit and the Luz uncertainty.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / ".claude/skills/whz-figure-style/scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from whz_plot import PALETTE, apply_style, fit_line, raw_trace, style_axes

from process_tierA_sop import USER_20261004_KNOTS, drift_corrected_nm
from reprocess_26_core_sop import (
    HANDLE_BP, MIN_DELTA_NM_PER_BP, MIN_INVERSION_FORCE, delta_nm_per_bp, x_ds,
)
from plot_26_core_sop import model_curve
from plot_whz_59 import PAGE_SIZE, draw_fec, draw_hist, draw_luz
from unzip_raw_points import load_raw_run
from whz_all_934 import process_file, raw_page
from tracking_quality import file_tracking

RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/20261004")
DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/"
            "Proecessed Traces/WHZ_20261004")

# Fields of view, from identical ROI coordinates; 051 and 052 share one field.
FOVS = {  # data ids per field; Fmax per file comes from the data
    "FOV01": [2, 3], "FOV02": [4, 5, 6, 7, 8], "FOV03": [9, 10, 11, 12, 13, 14],
    "FOV04": [15, 16, 17, 18, 19, 20], "FOV05": [21, 22, 23, 24, 25],
    "FOV06": [26, 27, 28, 29, 30, 31, 32], "FOV07": [33, 34, 35, 36, 37, 38, 39],
    "FOV08": [42], "FOV09": [44, 45, 46, 47, 48, 49, 50], "FOV10": [51, 52],
    "FOV11": [53, 54, 55, 56, 57, 58, 59],
}
FORCE_BIN_PN = 0.25        # 47 s ramps: 0.10 pN bins hold ~1 sample
# Force calibration for this batch (user, 2026-10-07): 1.7 / 1.0 / 0.7 mm =
# 18 / 30 / 46 pN. The SOP law (x1.15) gave 20.8 / 42.3 / 57.2 pN at those
# magnet positions; every force, anchor and Fmax label here uses the new one.
CALIBRATION = "user_20261004"

# Construct: YTrunk_2601_980bp (SnapGene, OneDrive DNA/DNA map from elab). The
# user normalises it to 980 bp; the SnapGene record is the 971 bp PCR product
# with BsaI sites at both ends. Two Widom 601 sites sit at 314-460 and 511-657.
# The construct is nearly symmetric (313 bp before the first site, 314 bp after
# the second), so the sites fall at the same Luz whichever end the fork enters.
CONSTRUCT_BP = 980.0
NUCLEOSOME_SITES_BP = ((314.0, 460.0), (511.0, 657.0))
# A bead counts as fully unzipped when its reference span is within this band
# of the dWLC span of a fully open 980 bp construct (797 nm, 14 -> 38 pN).
FULL_UNZIP_RATIO = (0.80, 1.25)
MIN_FULL_BEADS_PER_FOV = 3
F_OWN_STANDARD = 38.6      # the SOP's 38 pN full anchor plus its +-0.6 pN window
MIN_SPAN_NM = 200.0
# Upper bound on a physical span. With the SOP's dWLC parameters one construct
# goes from closed at 14 pN to fully open at 38 pN over ~1.1-1.2 um, and the
# longest trace the user accepted from 20260731 (D026_R063) spans 2.2 um.
# 2.5 um still admits a double tether. The first full run of this batch, with
# no upper bound, chose a 'reference' above 2.5 um for 169 of 477 beads (up to
# 207 um): beads lost or mistracked in their highest-force file, whose Z jumps
# made a garbage span that then squashed every other file of the bead.
MAX_SPAN_NM = 2500.0
# Reported, never used to reject. A lower-force file with a SMALLER own span
# than the reference is expected (it did not unzip fully); one with a span
# LARGER than the reference by this fraction unzipped further at lower force,
# which a single unchanged tether cannot do - it points to a tether change or
# a tracking artefact in one of the files.
REF_DISAGREE_FLAG = 0.30
LEVEL_JUMP_UM = 1.5                         # tracking_quality.TH["level_jump_um"]
ZERO_DONE_PN = 16.5                         # 14 pN zero anchor + its 2.5 pN reach
LUZ_GUARD_BP = (-300.0, 1.5 * CONSTRUCT_BP)
TOP_WINDOW_PN = 3.0
# The reference is the median of the reference file's cycles (user's choice,
# 2026-10-07). When those cycles differ by more than this fraction of their
# median, the median of two is the mean of two different states - one cycle
# stalled, or one cycle on a LUT level (D039_R026: 261 vs 1432 nm). Such a bead
# keeps its median scale but is flagged, does not count towards a full-unzip
# scale, and is left out of the stall statistics.
REF_CYCLE_SPREAD_FLAG = 0.30
CYC = {1: PALETTE["red"], 2: PALETTE["blue"]}


def bead_of(roi_key: str) -> str:
    return roi_key.split("_")[1]          # "D009_R012" -> "R012"


def expected_span_nm(zero_f: float = 14.0, full_f: float = 38.0) -> float:
    """dWLC span of a fully unzipped CONSTRUCT_BP construct between the anchors."""
    return float(HANDLE_BP * (x_ds(np.array([full_f])) - x_ds(np.array([zero_f])))[0]
                 + CONSTRUCT_BP * delta_nm_per_bp(np.array([full_f]))[0])


def normalise(ext_nm, force, z_zero, span, zero_f, full_f):
    """SOP section 4-5: 0-CONSTRUCT_BP normalisation and dWLC Luz."""
    norm = (np.asarray(ext_nm, float) - z_zero) * CONSTRUCT_BP / span
    f = np.asarray(force, float)
    dx = delta_nm_per_bp(f)
    ds = HANDLE_BP * (x_ds(f) - x_ds(np.full_like(f, zero_f)))
    e_full = (HANDLE_BP * (x_ds(np.array([full_f])) - x_ds(np.array([zero_f])))
              + CONSTRUCT_BP * delta_nm_per_bp(np.array([full_f])))[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        luz = (norm / (CONSTRUCT_BP / e_full) - ds) / dx
    bad = (f < MIN_INVERSION_FORCE) | (dx <= MIN_DELTA_NM_PER_BP) | ~np.isfinite(luz)
    return norm, np.where(bad, np.nan, luz)


def bead_references(anchors: pd.DataFrame) -> dict[str, dict]:
    """One length scale per bead, from its highest-force validly anchored file."""
    own = anchors[(anchors.fmax_pn >= F_OWN_STANDARD)
                  & anchors.raw_anchor_span_nm.between(MIN_SPAN_NM, MAX_SPAN_NM)
                  & (anchors.full_branch_fallback.fillna("") == "")
                  & anchors.anchors_tracked]
    refs = {}
    for bead, g in own.groupby("bead"):
        ref_file = int(g.loc[g.fmax_pn.idxmax(), "data_id"])
        rg = g[g.data_id == ref_file]
        span = float(np.median(rg.raw_anchor_span_nm))
        # comparable measurements only: the bead's other own-anchored files,
        # one median per file so a single noisy cycle cannot raise the flag
        lower = g[g.data_id != ref_file].groupby("data_id").raw_anchor_span_nm.median().values
        disagree = (float((lower.max() - span) / span)
                    if lower.size and span > 0 else np.nan)
        ratio = span / expected_span_nm(14.0, float(np.median(rg.full_force_pn)))
        kind = ("full" if FULL_UNZIP_RATIO[0] <= ratio <= FULL_UNZIP_RATIO[1]
                else "partial" if ratio < FULL_UNZIP_RATIO[0] else "over_long")
        spans = rg.raw_anchor_span_nm.values
        spread = float((spans.max() - spans.min()) / span) if span > 0 else np.nan
        refs[bead] = {"own_ref_span": span, "span": span, "data_id": ref_file,
                      "cycle_spread": spread,
                      "cycles_disagree": bool(spread > REF_CYCLE_SPREAD_FLAG),
                      "fmax_pn": float(rg.fmax_pn.iloc[0]),
                      "full_force_pn": float(np.median(rg.full_force_pn)),
                      "cycle_spans": [round(float(x)) for x in rg.raw_anchor_span_nm],
                      "disagree": disagree, "ratio": ratio, "kind": kind}
    return refs


def fov_full_scale(refs: dict) -> tuple[float, int]:
    """Median reference span of the beads in a field that unzipped fully."""
    full = [r["own_ref_span"] for r in refs.values()
            if r["kind"] == "full" and not r["cycles_disagree"]]
    return (float(np.median(full)), len(full)) if full else (np.nan, 0)


def rescale_incomplete(refs: dict, scale: float, n: int, scope: str) -> None:
    """A bead that never unzipped fully takes the field's full-unzip span.

    Its own span reaches only as far as it got - typically a nucleosome - so
    using it as the 980 bp scale would map a stall onto a complete opening.
    """
    for r in refs.values():
        if r["kind"] != "full":
            r["span"] = scale
            r["scale_note"] = (f"{r['kind'].replace('_', '-')} (own {r['own_ref_span']:.0f} nm "
                               f"= {r['ratio']:.2f}x); scale = {scope} full-unzip median "
                               f"{scale:.0f} nm, n={n}")
        else:
            r["scale_note"] = f"full unzip ({r['ratio']:.2f}x)"


def apply_anchor_rules(anc: pd.DataFrame, fmax: float, refs: dict) -> pd.DataFrame:
    a = anc.copy()
    a["own_span_nm"] = a.raw_anchor_span_nm
    for col in ("anchor_source", "ref_data_id", "ref_fmax_pn", "ref_disagree",
                "ref_kind", "ref_own_span_nm", "ref_ratio_to_construct",
                "ref_cycle_spread", "ref_cycles_disagree"):
        a[col] = "" if col in ("anchor_source", "ref_kind") else np.nan
    for i, r in a.iterrows():
        ref = refs.get(r.bead)
        zero_ok = np.isfinite(r.get("z_zero_raw_nm", np.nan))
        if r.tracking in ("no_bead", "lost_before_pull"):
            a.at[i, "anchor_valid"] = False
            a.at[i, "reason"] = f"bead {r.tracking.replace('_', ' ')} in this file"
            a.at[i, "anchor_source"] = "not normalised"
            continue
        if r.level_jump_bad:
            a.at[i, "anchor_valid"] = False
            a.at[i, "reason"] = f"LUT level jump {r.level_jump_um:+.1f} um in this file"
            a.at[i, "anchor_source"] = "not normalised"
            continue
        if r.tracking == "lost" and not r.t_loss_s > r.t_zero_done_s:
            a.at[i, "anchor_valid"] = False
            a.at[i, "reason"] = f"bead lost at {r.t_loss_s:.1f} s, before this cycle's zero"
            a.at[i, "anchor_source"] = "not normalised"
            continue
        if ref is None:
            a.at[i, "anchor_valid"] = False
            a.at[i, "reason"] = ("no file >= 38.6 pN with a valid anchor for this bead, "
                                 "so full unzipping is never confirmed")
            a.at[i, "anchor_source"] = "not normalised"
            continue
        if not zero_ok:
            a.at[i, "anchor_valid"] = False
            a.at[i, "reason"] = "no zero anchor in this file"
            a.at[i, "anchor_source"] = "not normalised"
            continue
        a.at[i, "raw_anchor_span_nm"] = round(ref["span"], 2)
        a.at[i, "full_force_pn"] = ref["full_force_pn"]
        a.at[i, "anchor_valid"] = ref["span"] >= MIN_SPAN_NM
        a.at[i, "reason"] = "" if ref["span"] >= MIN_SPAN_NM else (
            f"reference span {ref['span']:.0f} nm < {MIN_SPAN_NM:.0f} nm")
        a.at[i, "ref_kind"] = ref["kind"]
        a.at[i, "ref_own_span_nm"] = round(ref["own_ref_span"], 1)
        a.at[i, "ref_ratio_to_construct"] = round(ref["ratio"], 3)
        a.at[i, "ref_data_id"] = ref["data_id"]
        a.at[i, "ref_fmax_pn"] = round(ref["fmax_pn"], 1)
        a.at[i, "ref_disagree"] = ref["disagree"]
        a.at[i, "ref_cycle_spread"] = round(ref["cycle_spread"], 3)
        a.at[i, "ref_cycles_disagree"] = ref["cycles_disagree"]
        tag = "this file" if int(r.data_id) == ref["data_id"] else f"data_{ref['data_id']:03d}"
        a.at[i, "anchor_source"] = (f"ref {tag} (Fmax {ref['fmax_pn']:.0f} pN), "
                                    + ref.get("scale_note", "")
                                    + (f"; REF CYCLES DIFFER {ref['cycle_spread']:.0%}"
                                       if ref["cycles_disagree"] else ""))
    return a


def raw_table(run, anc: pd.DataFrame, fmax: float) -> pd.DataFrame:
    """Every raw sample of every ROI, normalised where the cycle's anchor holds.

    Samples after a lost bead's loss time are never normalised, and a cycle
    whose Luz at Fmax falls outside LUZ_GUARD_BP is dropped here; `anc` is
    updated in place so the anchor table and the pages say why.
    """
    pos = {(r.roi_key, int(r.cycle)): i for i, r in anc.iterrows()}
    anc["luz_at_fmax_bp"] = np.nan
    frames = []
    for roi in run.rois:
        key = f"D{run.data_id:03d}_R{roi}"
        ext = drift_corrected_nm(run.z_um[roi], run.time_s)
        for cyc, direction, seg in run.branches:
            i = pos.get((key, int(cyc)))
            r = anc.loc[i] if i is not None else None
            f = run.force_pn[seg]
            norm = luz = np.full(seg.size, np.nan)
            ok = bool(r is not None and r.anchor_valid and np.isfinite(r.z_zero_raw_nm))
            keep = np.full(seg.size, ok)
            if ok:
                norm, luz = normalise(ext[seg], f, r.z_zero_raw_nm, r.raw_anchor_span_nm,
                                      r.zero_force_pn, r.full_force_pn)
                if r.tracking == "lost":
                    keep &= run.time_s[seg] < r.t_loss_s
            frames.append(pd.DataFrame({
                "roi_key": key, "cycle": int(cyc), "branch": direction,
                "force_pn": np.round(f, 3), "extension_raw_nm": np.round(ext[seg], 1),
                "extension_normalized_nm": np.round(np.where(keep, norm, np.nan), 1),
                "luz_dwlc_bp": np.round(np.where(keep, luz, np.nan), 1),
                "normalized_valid": keep}))
    raw = pd.concat(frames, ignore_index=True)
    top = raw[raw.normalized_valid & (raw.force_pn >= fmax - TOP_WINDOW_PN)]
    luz_top = top.groupby(["roi_key", "cycle"]).luz_dwlc_bp.median()
    for (key, cyc), v in luz_top.items():
        i = pos[(key, int(cyc))]
        anc.at[i, "luz_at_fmax_bp"] = round(float(v), 1)
        hi = np.inf if anc.at[i, "ref_kind"] == "over_long" else LUZ_GUARD_BP[1]
        if np.isfinite(v) and LUZ_GUARD_BP[0] <= v <= hi:
            continue
        anc.at[i, "anchor_valid"] = False
        anc.at[i, "reason"] = (f"Luz at Fmax {v:.0f} bp outside {LUZ_GUARD_BP[0]:.0f}.."
                               f"{LUZ_GUARD_BP[1]:.0f}: tracking broken in this cycle")
        anc.at[i, "anchor_source"] = "not normalised"
        sel = (raw.roi_key == key) & (raw.cycle == cyc)
        raw.loc[sel, ["extension_normalized_nm", "luz_dwlc_bp"]] = np.nan
        raw.loc[sel, "normalized_valid"] = False
    return raw


def file_pages(pdf_path, data_id, fmax, raw, anc):
    index = []
    with PdfPages(pdf_path) as pdf:
        for roi in sorted(raw.roi_key.unique()):
            p = raw[raw.roi_key == roi]
            a = anc[anc.roi_key == roi]
            fig, axes = plt.subplots(1, 3, figsize=PAGE_SIZE)
            valid = a[a.anchor_valid]
            if len(valid) and p.normalized_valid.any():
                pv = p[p.normalized_valid]
                draw_fec(axes[0], pv, valid, full_bp=CONSTRUCT_BP)
                draw_luz(axes[1], pv, full_bp=CONSTRUCT_BP, bands=NUCLEOSOME_SITES_BP)
                draw_hist(axes[2], pv, full_bp=CONSTRUCT_BP, bands=NUCLEOSOME_SITES_BP)
                note = "; ".join(f"C{int(r.cycle)}: {r.anchor_source}" for r in valid.itertuples())
                dropped = a[~a.anchor_valid]
                if len(dropped):
                    note += "; " + "; ".join(f"C{int(r.cycle)} excluded: {r.reason}"
                                             for r in dropped.itertuples())
            else:
                raw_page(fig, axes, p, a)
                note = "; ".join(f"C{int(r.cycle)} excluded: {r.reason}" for r in a.itertuples())
            trk = a.tracking.iloc[0] if len(a) else "unknown"
            if trk == "lost":
                trk += f" at {a.t_loss_s.iloc[0]:.1f} s"
            fig.suptitle(f"20261004_{roi}   Fmax {fmax:.1f} pN · tracking: {trk} · "
                         f"{len(p)} raw points · {note}", fontsize=6.0, y=1.03)
            fig.tight_layout(w_pad=1.6)
            pdf.savefig(fig, dpi=300, bbox_inches="tight")
            plt.close(fig)
            index.append({"pdf_page": len(index) + 1, "roi_key": roi,
                          "n_valid_cycles": int(a.anchor_valid.sum()), "n_cycles": len(a)})
    return pd.DataFrame(index)


def titration_pages(pdf_path, fov, files, raws, ancs, fmaxs):
    """One page per bead: its FEC at every maximum force, same axes throughout."""
    order = sorted(files, key=lambda d: fmaxs[d])
    beads = sorted({bead_of(k) for d in files for k in raws[d].roi_key.unique()})
    ncol = 4
    nrow = int(np.ceil(len(order) / ncol))
    xmax = max(fmaxs.values()) + 1
    with PdfPages(pdf_path) as pdf:
        for bead in beads:
            # limits are set per panel rather than shared, so every visible panel
            # keeps its own tick labels when the bottom row is only part-filled
            fig, axes = plt.subplots(nrow, ncol, figsize=(7.09, 2.0 * nrow + 0.3),
                                     squeeze=False)
            for ax in axes.flat:
                ax.set_visible(False)
            for k, d in enumerate(order):
                ax = axes.flat[k]
                ax.set_visible(True)
                key = f"D{d:03d}_{bead}"
                p = raws[d][raws[d].roi_key == key]
                a = ancs[d][(ancs[d].roi_key == key) & ancs[d].anchor_valid]
                norm = p.normalized_valid.any()
                for cyc in (2, 1):
                    for br, filled in (("down", False), ("up", True)):
                        g = p[(p.cycle == cyc) & (p.branch == br)]
                        if norm:
                            g = g[g.normalized_valid]
                        if g.empty:
                            continue
                        raw_trace(ax, g.force_pn,
                                  g.extension_normalized_nm if norm else g.extension_raw_nm,
                                  color=CYC.get(cyc, PALETTE["dark_grey"]), filled=filled,
                                  marker_area=1.2, alpha=0.5, edge_width=0.25, rasterized=True)
                if norm and len(a):
                    r = a.iloc[0]
                    grid = np.linspace(5.0, fmaxs[d], 200)
                    for luz, ls in ((0.0, "--"), (CONSTRUCT_BP, "-")):
                        fit_line(ax, grid, model_curve(grid, luz, r.zero_force_pn,
                                                       r.full_force_pn, CONSTRUCT_BP),
                                 color=PALETTE["light_grey"], linestyle=ls, linewidth=0.7,
                                 zorder=1)
                ax.set_xlim(0, xmax)
                ax.set_ylim(-300, 1.3 * CONSTRUCT_BP)
                trk = ancs[d][ancs[d].roi_key == key].tracking
                trk = trk.iloc[0] if len(trk) else "unknown"
                tag = ("" if norm else " (raw nm)") + ("" if trk == "tracked" else f" [{trk}]")
                is_ref = bool(len(a)) and float(a.ref_data_id.iloc[0]) == d
                ax.set_title(f"data_{d:03d}  Fmax {fmaxs[d]:.0f} pN{tag}"
                             + ("  [length ref]" if is_ref else ""), fontsize=6.5)
                style_axes(ax, xlabel="Force (pN)",
                           ylabel="Norm. ext. (nm)" if k % ncol == 0 else None)
                if k % ncol:
                    ax.tick_params(labelleft=False)
                ax.tick_params(labelsize=6.5)
                ax.xaxis.label.set_size(7)
                ax.yaxis.label.set_size(7)
            fig.suptitle(f"20261004 {fov} bead {bead}: titration across {len(order)} files",
                         fontsize=7.5, y=1.0)
            fig.tight_layout(h_pad=0.6, w_pad=0.4)
            pdf.savefig(fig, dpi=200, bbox_inches="tight")
            plt.close(fig)
    return len(beads)


def main() -> int:
    apply_style()
    DEST.mkdir(parents=True, exist_ok=True)

    # pass 1: anchors, tracking and a reference for every bead of every field
    binned, ancs, fmaxs, refs_by_fov = {}, {}, {}, {}
    for fov, files in FOVS.items():
        for d in files:
            run = load_raw_run(RAW, d, CALIBRATION)
            fmaxs[d] = float(run.force_pn.max())
            # When each cycle has finished with its anchors: the down branch has
            # left the full-anchor window (38 pN, up to 4 pN reach). A bead lost
            # after that point still has a usable anchor for that cycle; on these
            # ramps most "losses" happen in the final near-zero-force dwell.
            anchors_done, zero_done = {}, {}
            for cyc, branch, seg in run.branches:
                if branch == "down":
                    below = np.flatnonzero(run.force_pn[seg] <= 33.5)
                    anchors_done[int(cyc)] = float(run.time_s[seg][below[0]] if below.size
                                                   else run.time_s[seg][-1])
                else:
                    above = np.flatnonzero(run.force_pn[seg] >= ZERO_DONE_PN)
                    zero_done[int(cyc)] = float(run.time_s[seg][above[0]] if above.size
                                                else run.time_s[seg][-1])
            pts, anc = process_file(d, RAW, FORCE_BIN_PN, CALIBRATION)
            trk = file_tracking(RAW / f"data_{d:03d}.h5")
            roi_of = anc.roi_key.str.split("_R").str[1]
            for col, key in (("tracking", "status"), ("tracking_bad_fraction", "bad_fraction"),
                             ("t_loss_s", "t_loss_s"), ("level_jump_um", "level_jump_um")):
                default = "unknown" if key == "status" else np.nan
                anc[col] = roi_of.map(lambda r, k=key, dflt=default: trk.get(r, {}).get(k, dflt))
            anc["t_anchors_done_s"] = anc.cycle.map(lambda c: anchors_done.get(int(c), np.nan))
            anc["t_zero_done_s"] = anc.cycle.map(lambda c: zero_done.get(int(c), np.nan))
            # a lost bead's level jump is measured before the loss; it is the same
            # LUT artefact as in a bead that stays tracked
            anc["level_jump_bad"] = anc.level_jump_um.abs() > LEVEL_JUMP_UM
            anc["anchors_tracked"] = ((anc.tracking == "tracked") | (
                (anc.tracking == "lost") & (anc.t_loss_s > anc.t_anchors_done_s))
            ) & ~anc.level_jump_bad
            anc["data_id"], anc["fov"], anc["fmax_pn"] = d, fov, fmaxs[d]
            anc["bead"] = anc.roi_key.map(bead_of)
            binned[d], ancs[d] = pts, anc
        refs_by_fov[fov] = bead_references(pd.concat([ancs[d] for d in files],
                                                     ignore_index=True))
        print(f"{fov}: anchors done for {len(files)} files", flush=True)

    # length scale for beads that never unzipped fully: field median, else batch
    batch_full = [r["own_ref_span"] for refs in refs_by_fov.values()
                  for r in refs.values() if r["kind"] == "full" and not r["cycles_disagree"]]
    batch_scale = float(np.median(batch_full)) if batch_full else np.nan
    scale_log = {}
    for fov, refs in refs_by_fov.items():
        scale, n = fov_full_scale(refs)
        scope = fov
        if n < MIN_FULL_BEADS_PER_FOV:
            scale, n, scope = batch_scale, len(batch_full), "batch"
        rescale_incomplete(refs, scale, n, scope)
        kinds = pd.Series([r["kind"] for r in refs.values()]).value_counts().to_dict()
        scale_log[fov] = {"scale_nm": round(scale, 1) if np.isfinite(scale) else None,
                          "scope": scope, "n_full": n, "beads_by_kind": kinds}
        print(f"{fov}: references {kinds}; incomplete beads scaled by {scope} "
              f"full-unzip median {scale:.0f} nm (n={n})", flush=True)

    # pass 2: normalise, draw, write
    summary, all_anc = [], []
    for fov, files in FOVS.items():
        refs = refs_by_fov[fov]
        raws = {}
        for d in files:
            run = load_raw_run(RAW, d, CALIBRATION)
            ancs[d] = apply_anchor_rules(ancs[d], fmaxs[d], refs)
            raws[d] = raw_table(run, ancs[d], fmaxs[d])
            raws[d].to_csv(DEST / f"raw_points_data_{d:03d}.csv.gz", index=False)
            binned[d].to_csv(DEST / f"binned_points_data_{d:03d}.csv.gz", index=False)
            n_roi = raws[d].roi_key.nunique()
            idx = file_pages(DEST / f"WHZ_20261004_data_{d:03d}_{n_roi}_ROIs.pdf",
                             d, fmaxs[d], raws[d], ancs[d])
            idx.to_csv(DEST / f"WHZ_20261004_data_{d:03d}_page_index.csv", index=False)
            a = ancs[d]
            tr = a.drop_duplicates("roi_key").tracking.value_counts()
            kinds = a[a.anchor_valid].drop_duplicates("roi_key").ref_kind.value_counts()
            nv = int(a.groupby("roi_key").anchor_valid.any().sum())
            summary.append({"fov": fov, "data_id": d, "fmax_pn": round(fmaxs[d], 1),
                            "n_roi": n_roi, "roi_with_valid_cycle": nv,
                            "valid_cycles": int(a.anchor_valid.sum()), "n_cycles": len(a),
                            "tracked": int(tr.get("tracked", 0)), "lost": int(tr.get("lost", 0)),
                            "absent": int(tr.get("no_bead", 0) + tr.get("lost_before_pull", 0)),
                            "level_jump": int(tr.get("level_jump_artefact", 0)),
                            "beads_full_unzip": int(kinds.get("full", 0)),
                            "beads_partial_unzip": int(kinds.get("partial", 0)),
                            "beads_over_long": int(kinds.get("over_long", 0))})
            all_anc.append(a)
            print(f"{fov} data_{d:03d} Fmax {fmaxs[d]:5.1f} pN: {n_roi} ROI -> {nv} normalised; "
                  f"tracked {int(tr.get('tracked', 0))}", flush=True)
        if len(files) > 1:
            nb = titration_pages(DEST / f"WHZ_20261004_{fov}_titration_by_bead.pdf",
                                 fov, files, raws, ancs, fmaxs)
            print(f"{fov}: titration PDF, {nb} beads x {len(files)} files", flush=True)

    s = pd.DataFrame(summary)
    s.to_csv(DEST / "WHZ_20261004_summary.csv", index=False)
    anchors = pd.concat(all_anc, ignore_index=True)
    anchors.to_csv(DEST / "per_cycle_anchors.csv", index=False)
    beads = anchors[anchors.anchor_valid].drop_duplicates(["fov", "bead"])
    print(f"\nTOTAL {s.n_roi.sum()} ROI-file instances, {s.roi_with_valid_cycle.sum()} "
          f"normalised; beads with a scale {len(beads)}/{anchors.drop_duplicates(['fov','bead']).shape[0]}: "
          f"{beads.ref_kind.value_counts().to_dict()}; batch full-unzip median {batch_scale:.0f} nm "
          f"(n={len(batch_full)}), dWLC expectation {expected_span_nm():.0f} nm")
    (DEST / "manifest.json").write_text(json.dumps({
        "sop": "20260811_unzipnucleosomedata_complete_processing_pipeline.docx",
        "raw_dir": str(RAW), "fields_of_view": FOVS,
        "construct": {"name": "YTrunk_2601_980bp", "normalised_to_bp": CONSTRUCT_BP,
                      "snapgene_length_bp": 971,
                      "nucleosome_601_sites_bp": NUCLEOSOME_SITES_BP,
                      "expected_full_span_nm_14_to_38pN": round(expected_span_nm(), 1)},
        "force_calibration": {"name": CALIBRATION, "knots_shift_mm_force_pn": USER_20261004_KNOTS,
                              "interpolation": "piecewise log-linear; end segments extended "
                                               "outside 0.7-1.7 mm",
                              "replaces": "SOP 20260724 law x1.15 (20.8/42.3/57.2 pN at "
                                          "1.7/1.0/0.7 mm)"},
        "force_bin_pn": FORCE_BIN_PN,
        "anchor_rules": {
            "anchors": "estimated on raw samples; zero >=5 samples within 2.5 pN of 14 pN, "
                       "full >=3 samples within 4 pN of 38 pN; Theil-Sen to the anchor force",
            "reference": "per bead, its highest-Fmax file (>=38.6 pN) in which it is tracked "
                         "and has a 200-2500 nm span; median of that file's valid cycles",
            "full_unzip_ratio": FULL_UNZIP_RATIO,
            "incomplete_beads": "partial or over-long references take the field's full-unzip "
                                f"median span; fewer than {MIN_FULL_BEADS_PER_FOV} full beads in "
                                "the field -> batch median",
            "length_scale_applied_to": "every file of the bead; zero from each file and cycle",
            "no_reference": "not normalised; drift-corrected extension shown",
            "never_normalised": "files where the bead is absent or lost before the pull",
            "reference_requires_tracking": "per cycle: tracked for the whole file, or lost only "
                                           "after that cycle's down branch left the 38 pN "
                                           "anchor window (tracking_quality.py, 20260913 criteria)",
            "ref_disagree_flag": REF_DISAGREE_FLAG,
            "tracking_gate_every_file": "level jump > 1.5 um -> file not normalised; lost -> "
                                        "samples after the loss not normalised, cycle dropped "
                                        "if lost before its zero window (16.5 pN up)",
            "luz_guard_bp": LUZ_GUARD_BP,
            "ref_cycle_spread_flag": REF_CYCLE_SPREAD_FLAG,
            "ref_cycles_disagree": "reference cycles differ by more than the flag: median kept "
                                   "(user's rule), bead flagged, not used for a full-unzip "
                                   "scale",
            "luz_guard": "cycle dropped if median Luz within 3 pN of Fmax is outside the "
                         "guard (upper bound not applied to over-long beads)"},
        "length_scales_by_fov": scale_log,
        "batch_full_unzip_median_nm": round(batch_scale, 1),
        "display": "raw samples; bins feed the anchor fit and Luz uncertainty only",
        "style": "whz-figure-style; cycle 1 red #E64B35, cycle 2 blue #3C5488; "
                 "601 sites shaded on the Luz panels",
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
