# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: .venv
#     language: python
#     name: python3
# ---

# %%
from pathlib import Path
import matplotlib.pyplot as plt
from trace_editor.core.h5_data import H5Data
import numpy as np

# %matplotlib inline

# %%
# --- configure here ---
date = 20260706
FOLDER      = Path("/Users/barton/My Drive/data") / str(date)

FILE_RANGE  = range(20, 27)                        # data_001.h5 to data_018.h5
FILES       = [FOLDER / f"data_{i:03d}.h5" for i in FILE_RANGE]
ROI_RANGE   = [10] #range(53, 54)                         # e.g. range(0,10) or [2,5,7]
X_CH        = "Focus (mm)" 
X_CH      = "Time (s)"
# X_CH      = "Intensity"
# X_CH        = "Shift (mm)"
Y_CH        = "Z (um)"
Idep_var    = "Shift (mm)"
LEGEND_KEY  =  "Motors/Focus_mm" # #"Imaging/exposure_time_ms" 'Motors/Shift_mm'
# ----------------------

def get_setting(path, group_key):
    """Read a Group/Key value from the Settings table, e.g. 'Motors/Shift_mm'."""
    import h5py
    group, key = group_key.split("/", 1)
    with h5py.File(path, "r") as f:
        for row in f["Settings"][:]:
            if row[0].decode() == group and row[1].decode() == key:
                return row[2].decode()
    return "?"

print(f"{len(FILES)} files")
for f in FILES:
    print(" ", f.name)

# %%
fig, ax = plt.subplots(figsize=(10, 4))
ax_idv = ax.twinx()

cmap = plt.get_cmap('plasma')
label_values = np.array([float(get_setting(path, LEGEND_KEY)) for path in FILES])

for path in FILES:
    h5 = H5Data(path)
    all_rois = h5.list_rois()
    focus_val = float(get_setting(path, LEGEND_KEY))*1000
    label_val = f"{float(get_setting(path, LEGEND_KEY))*1000}"

    for roi_idx in ROI_RANGE:
        roi = f"{roi_idx:03d}"
        if roi not in all_rois:
            continue

        x = h5.get_channel(X_CH, roi) #* 1000
        y = h5.get_channel(Y_CH, roi)
        IDV = h5.get_channel(Idep_var, roi) #* 1000

        lbl = f"{label_val} roi{roi}"
        # y += label_val * 1000
        ax.scatter(x, y - y[0], label=lbl, s=2)
        # ax_idv.plot(x, -IDV, color='grey', alpha=0.7)
# ax.set_ylim(0, 30)
plt.title(f"{date}, fileN={FILE_RANGE}, ROI_N={ROI_RANGE}, Key={LEGEND_KEY}")

ax.set_xlabel(X_CH)
ax.set_ylabel(Y_CH)
# ax_idv.set_ylabel(Idep_var)
ax.set_ylim(-2)
ax.legend(fontsize=9, ncol=2, loc=1, title=f'{LEGEND_KEY}')

plt.tight_layout()
plt.show()
