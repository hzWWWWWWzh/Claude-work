"""Helper for merge_selected_20261013.py: run luzB_peakmap_by_force_20261009 unchanged on the 20261013 selection."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import process_20261009 as p
p.RAW = Path("/Users/whz/Desktop/DATA MT/unzipping/20261013")
p.DEST = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces/WHZ_20261013")
import ladder_norm1028_20261009 as L
L.SEL_CSV = p.DEST / "selection_101_tethers.csv"
import matplotlib.figure as mf
_orig = mf.Figure.suptitle
def _suptitle(self, t, *a, **k):          # the module's title text names the 1009 batch
    return _orig(self, t.replace("20261009, your 96 selected tethers", "20261013 (+1 mM Mg2+), your 101 selected tethers"), *a, **k)
mf.Figure.suptitle = _suptitle
import luzB_peakmap_by_force_20261009 as M
raise SystemExit(M.main())
