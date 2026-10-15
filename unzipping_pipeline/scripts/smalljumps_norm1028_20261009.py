#!/usr/bin/env python3
"""20261009: the small-jump analysis (smalljumps_modelB_20261009.py) re-expressed on the
user's 0-1028 scale (user, 2026-10-11: "使用我的方法吧").

Every detected step keeps its time and size in nm; its size on the 0-1028 scale is
1028 * step_nm / span (span = the tether's 1028 point - zero, norm1028_selected_96_peakmap.csv).
BIG = forward >= 100, SMALL = |step| < 100 on that scale. Open-fork plateaus >= 18 pN of the
user's 96 tethers; background = the same tethers' 10 pN plateaus (closed fork, nothing
can open). Near a big jump = within 1 s after/before it (plateau arrival counted as a big
event as in the model-B analysis); far = > 2 s from every big event.
Detection: 80 % recall needs ~60 nm (verify: A80 ~54 bp B at ~1.1 nm/bp), i.e.
1028 * 60 / span on this scale - for spans < 617 nm every 'small' (< 100) jump is below it.
Output: smalljumps_modelB/norm1028_scale/ (tables + one figure).
"""
from pathlib import Path
import numpy as np
import pandas as pd

D = Path("/Users/whz/Desktop/PHD project/Unzipping/Unzipping nucleoome/Proecessed Traces/WHZ_20261009")
IN = D / "smalljumps_modelB"
OUT = IN / "norm1028_scale"
OUT.mkdir(exist_ok=True)
A80_NM = 60.0

span = pd.read_csv(D / "norm1028_selected_96_peakmap.csv")
span["tether"] = span.data_id.astype(str) + "_" + span.bead
span = span.set_index("tether").span_nm
s = pd.read_csv(IN / "steps_all_flagged.csv.gz")
s = s[(s.set == "user96") & s.use & ~s.copy_suspect & s.tether.isin(span.index)].copy()
s["span_nm"] = s.tether.map(span)
s = s[s.span_nm > 0]
s["step_n"] = 1028 * s.step_nm / s.span_nm
s["cls_n"] = np.where(s.step_n >= 100, "big_fwd", np.where(s.step_n <= -100, "big_back",
                      np.where(s.step_n >= 0, "small_fwd", "small_back")))
s["detectable"] = s.step_nm.abs() >= A80_NM
pl = pd.read_csv(IN / "plateau_summary.csv")
pl = pl[pl.set == "user96"].copy()
pl["tether"] = pl.data_id.astype(str) + "_" + pl.bead
pl = pl[pl.tether.isin(span.index)]
qc = pd.read_csv(IN / "plateaus_qc.csv")[["data_id", "bead", "step", "hold_s", "passes"]]
pl = pl.merge(qc, on=["data_id", "bead", "step"], how="left")
key = ["data_id", "bead", "step"]

rows, near_rows = [], []
for f, g in pl.groupby("force_pn"):
    open_ = g[g.opened] if f >= 18 else g
    expo = open_.exp_open_s.sum() if f >= 18 else open_.exp_closed_s.sum()
    st = s.merge(open_[key], on=key)
    st = st[st.open_fork] if f >= 18 else st[~st.open_fork]
    sm = st[st.cls_n.str.startswith("small")]
    per_pl = sm.groupby(key).size()
    P = len(per_pl) / max(len(open_), 1)
    smd = sm[sm.detectable]
    rows.append({"force_pn": f, "fork": "open" if f >= 18 else "closed (background)",
                 "n_plateaus": len(open_), "exposure_s": round(expo, 1),
                 "n_small": len(sm), "n_small_fwd": int((sm.cls_n == "small_fwd").sum()),
                 "n_big_fwd": int((st.cls_n == "big_fwd").sum()),
                 "P_plateau_ge1_small": round(P, 3), "rate_small_per_s": round(len(sm) / expo, 4) if expo else np.nan,
                 "n_small_detectable(>=60nm)": len(smd),
                 "rate_small_detectable_per_s": round(len(smd) / expo, 4) if expo else np.nan})
tab = pd.DataFrame(rows)
# near / far around big forward jumps (new classification) and the arrival, open fork >= 18 pN
o = s[(s.force_pn >= 18) & s.open_fork].copy()
big = o[o.cls_n == "big_fwd"][key + ["t_s"]]
small = o[o.cls_n.str.startswith("small")]
hold = pl.set_index(key).hold_s.to_dict()
near_n = after_n = before_n = far_n = 0
near_t = far_t = 0.0
for k_, g in small.groupby(key):
    pass
for k_, p_ in pl[(pl.force_pn >= 18) & pl.opened].set_index(key).iterrows():
    h = hold.get(k_, np.nan)
    if not np.isfinite(h):
        continue
    tb = np.r_[0.3, big[(big.data_id == k_[0]) & (big.bead == k_[1]) & (big.step == k_[2])].t_s.values]
    grid = np.arange(0.3, h, 0.01)
    d = np.min(np.abs(grid[:, None] - tb[None, :]), axis=1)
    near_t += (d <= 1.0).sum() * 0.01
    far_t += (d > 2.0).sum() * 0.01
    ts = small[(small.data_id == k_[0]) & (small.bead == k_[1]) & (small.step == k_[2])].t_s.values
    if ts.size:
        ds = np.min(np.abs(ts[:, None] - tb[None, :]), axis=1)
        near_n += int((ds <= 1.0).sum())
        far_n += int((ds > 2.0).sum())
enr = pd.DataFrame([{"window": "within 1 s of a big jump or the arrival", "n_small": near_n,
                     "exposure_s": round(near_t, 1), "rate_per_s": round(near_n / near_t, 4)},
                    {"window": "> 2 s from every big event", "n_small": far_n,
                     "exposure_s": round(far_t, 1), "rate_per_s": round(far_n / far_t, 4)}])
enr["ratio_to_far"] = (enr.rate_per_s / enr.rate_per_s.iloc[1]).round(2)
tab.to_csv(OUT / "small_jump_probability_by_force_norm1028.csv", index=False)
enr.to_csv(OUT / "enrichment_near_big_norm1028.csv", index=False)
thr = (1028 * A80_NM / span[span > 0])
print(tab.to_string(index=False))
print(enr.to_string(index=False))
print("detection limit on the 0-1028 scale (80 % recall, 60 nm): per-tether quartiles",
      np.percentile(thr, [25, 50, 75]).round(0), "; tethers where it is >= 100:", int((thr >= 100).sum()), "of", thr.size)
amp = s[(s.force_pn >= 18) & s.open_fork & s.cls_n.str.startswith("small") & s.detectable].step_n
print("detectable small-jump sizes (0-1028 units), quartiles:", np.percentile(amp.abs(), [25, 50, 75]).round(0), "n", amp.size)
