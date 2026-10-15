#!/usr/bin/env python3
"""20261013 (Mg2+ added): every ROI of every titration file, one page each, for the user to pick.

User (2026-10-13): "处理1013的数据。我自己来挑选" / "今天加了mg2+，测的不是很好".
Same protocol as 20261009 (10 -> 39.9 pN ladder, 10 s holds, 15 s rests at 1.6 pN; some files
end with a 5 s 0.70 mm / ~46 pN step), so the 20261009 pipeline is reused unchanged, only pointed at
the 20261013 folder (process_20261009.RAW / DEST are reassigned here; no 1009 script is edited):
  * zero = maximum extension of the 10 pN step (1 s running median after a Hampel filter);
  * drift = line through the closed-fork rest dwells;
  * 1028 = furthest sustained extension on the highest step the bead SURVIVED (a step after a loss,
    fly-away, tracking failure or an unclosed rest is not used) - ladder_norm1028_20261009.bead_scale;
  * page = ladder_norm1028_20261009.draw with SCALE=norm1028: raw trace + FEC, the 0-1028 trace + its
    peak map on the same axis, force; the flags under the page (lost, skipped steps, small / over-long span,
    LUT level jumps) are those of the 20261009 selection pages.
No ROI is pre-selected or dropped; delete the pages you do not want, as for 20261009.

Outputs (Proecessed Traces/WHZ_20261013/select_pages/):
  WHZ_20261013_data_NNN_<n>_ROIs.pdf     one page per ROI, in ROI order
  select_page_index.csv                  pdf, page, data_id, bead, tracking, zero, 1028 step / force, span, flags
  files_skipped.csv                      unreadable (0-byte) files and non-titration files

    .venv/bin/python select_pages_20261013.py [data ids ...]
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MODEL", "B")
os.environ.setdefault("SCALE", "norm1028")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import process_20261009 as p  # noqa: E402

p.RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/20261013")
p.DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces/WHZ_20261013")

import ladder_norm1028_20261009 as L  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from tracking_quality import file_tracking  # noqa: E402

BATCH = "20261013"
OUT = p.DEST / "select_pages"
LEVELS = sorted({10.0, 18.0, 25.9, 26.9, 27.9, 30.0, 32.2, 34.6, 37.1, 39.9, 46.0})


def classify() -> tuple[list[int], pd.DataFrame]:
    """Titration files (with a step script); unreadable or other files are listed, not used."""
    titr, other = [], []
    for f_ in sorted(p.RAW.glob("data_*.h5")):
        d = int(f_.stem.split("_")[1])
        if f_.stat().st_size == 0:
            other.append({"data_id": d, "kind": "0 bytes (copy incomplete?)"})
            continue
        try:
            with h5py.File(f_) as h:
                steps = "Steps" in h
                tr = h["Traces"]
                n = len(tr["Z (um)"]) if "Z (um)" in tr else 0
                dur = float(tr["Time (s)"][-1] - tr["Time (s)"][0])
        except OSError as e:
            other.append({"data_id": d, "kind": f"unreadable: {e}"[:80]})
            continue
        if steps:
            titr.append(d)
        else:
            lut = f_.with_suffix(".lut").exists()
            other.append({"data_id": d, "kind": "LUT" if lut else "no step script (low-force check)",
                          "n_rois": n, "duration_s": round(dur, 1)})
    return titr, pd.DataFrame(other)


def fallback_page(run, ext, rec, title):
    """No 10 pN zero (bead lost or untracked before the end of the first step): raw trace and force
    only, so the page is still there to be kept or deleted."""
    from whz_plot import PALETTE, style_axes
    t = run["time_s"] - run["time_s"][0]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.09, 3.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    a1.scatter(t, ext, s=0.5, lw=0, color=PALETTE["dark_grey"], rasterized=True)
    a2.scatter(t, run["force_pn"], s=0.5, lw=0, color=PALETTE["dark_grey"], rasterized=True)
    style_axes(a1, ylabel="Extension (nm)")
    style_axes(a2, xlabel="Time (s)", ylabel="Force (pN)")
    fig.suptitle(title + "   - no 10 pN zero, not normalised", fontsize=7.5)
    fig.text(0.01, -0.02, "! " + " | ".join(rec["flags"] or ["no zero"]), fontsize=6.5)
    return fig


def main() -> int:
    p.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    titr, other = classify()
    other.to_csv(OUT / "files_skipped.csv", index=False)
    want = [int(a) for a in sys.argv[1:]] or titr
    colours = p.level_colours(LEVELS)
    old = pd.read_csv(OUT / "select_page_index.csv") if (OUT / "select_page_index.csv").exists() else pd.DataFrame()
    rows = []
    for d in [x for x in titr if x in want]:
        run = p.load(d)
        steps = p.segment_steps(run)
        trk = file_tracking(p.RAW / f"data_{d:03d}.h5")
        pdf_path = OUT / f"WHZ_{BATCH}_data_{d:03d}_{len(run['rois'])}_ROIs.pdf"
        with PdfPages(pdf_path) as pdf:
            for k, roi in enumerate(run["rois"], start=1):
                tq = trk.get(roi, {})
                ext, rec = L.bead_scale(run, steps, roi, tq)
                if np.isfinite(rec["end_nm"]):
                    rec["end_luz"] = float(L.model_luz(np.array([rec["end_nm"] - rec["zero_nm"]]),
                                                       np.array([rec["ref_force_pn"]]))[0])
                title = (f"{BATCH} (+Mg2+) data_{d:03d} R{roi}   page {k}/{len(run['rois'])}   "
                         f"tracking {tq.get('status', '?')}")
                if np.isfinite(rec["zero_nm"]):
                    fig, peaks = L.draw(run, steps, roi, ext, rec, colours, title)
                else:
                    fig, peaks = fallback_page(run, ext, rec, title), []
                pdf.savefig(fig, dpi=200, bbox_inches="tight")
                plt.close(fig)
                rows.append({"pdf": pdf_path.name, "page": k, "data_id": d, "bead": f"R{roi}",
                             "tracking": tq.get("status", "?"),
                             "t_loss_s": round(rec["t_loss"] - run["time_s"][0], 1) if np.isfinite(rec["t_loss"]) else np.nan,
                             "zero_10pN_max_nm": round(rec["zero_nm"], 1), "ref_step": rec["ref_step"],
                             "ref_force_pn": rec["ref_force_pn"], "end_max_nm": round(rec["end_nm"], 1),
                             "span_nm": round(rec["span_nm"], 1), "end_luzB_bp": round(rec.get("end_luz", np.nan), 1),
                             "skipped_steps": rec["skipped"], "max_held_s": round(rec.get("max_held_s", np.nan), 2),
                             "peaks_norm1028": "; ".join(f"{x:.0f}" for x in peaks), "drift": rec["drift"],
                             "flags": " | ".join(rec["flags"])})
        print(f"data_{d:03d}: {len(run['rois'])} pages -> {pdf_path.name}", flush=True)
    new = pd.DataFrame(rows)
    if len(old):
        old = old[~old.data_id.isin(new.data_id.unique())]
    idx = pd.concat([old, new], ignore_index=True).sort_values(["data_id", "page"])
    idx.to_csv(OUT / "select_page_index.csv", index=False)
    print(f"{len(new)} pages in {new.pdf.nunique() if len(new) else 0} PDFs; skipped files: "
          f"{other.kind.value_counts().to_dict() if len(other) else {}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
