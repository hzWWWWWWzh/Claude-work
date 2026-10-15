#!/usr/bin/env python3
"""20261013 (+1 mM Mg2+): merge the user's 101 selected traces with the ORIGINAL method.

User (2026-10-14): "对于整合和algin所有的trace，你还是别过多处理了。就按照原来的方法，帮我融合所有我选择的trace就可以了".
No re-anchoring: zero = max extension of the 10 pN step, 1028 = furthest sustained extension on the highest
step the bead survived (ladder_norm1028_20261009.bead_scale), exactly as on the selection pages.
  1) WHZ_20261013_selected_101_traces.pdf - the user's kept pages, unchanged, in their order, with bookmarks
     (from _user_selection_backup_101/);
  2) per-force pooled peak maps of the 101 tethers with luzB_peakmap_by_force_20261009 (unchanged) on the
     user's 0-1028 scale (norm1028_by_force/) and on Luz model B (luz_modelB/), process_20261009 pointed at 20261013.
    .venv/bin/python merge_selected_20261013.py
"""
import os, sys, subprocess
from pathlib import Path
import fitz, pandas as pd

HERE = Path(__file__).resolve().parent
D = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces/WHZ_20261013")
sel = pd.read_csv(D / "selection_101_tethers.csv")

out = fitz.open()
toc = []
for pdf, g in sel.groupby("pdf", sort=True):
    src = fitz.open(D / "_user_selection_backup_101" / pdf)
    toc.append([1, f"data_{int(g.data_id.iloc[0]):03d} ({len(g)})", out.page_count + 1])
    for r in g.sort_values("page_in_trimmed").itertuples():
        toc.append([2, f"{r.bead}  span {r.span_nm:.0f} nm", out.page_count + 1])
        out.insert_pdf(src, from_page=int(r.page_in_trimmed) - 1, to_page=int(r.page_in_trimmed) - 1)
out.set_toc(toc)
merged = D / f"WHZ_20261013_selected_{len(sel)}_traces.pdf"
out.save(merged, garbage=4, deflate=True)
print(f"{out.page_count} pages -> {merged}")
sel.to_csv(D / f"norm1028_selected_{len(sel)}.csv", index=False)

if "--no-peakmaps" not in sys.argv:
    for scale in ("norm1028", "luzB"):
        env = dict(os.environ, SCALE=scale, MODEL="B")
        subprocess.run([sys.executable, str(HERE / "_peakmap_1013_runner.py")], env=env, check=True)
