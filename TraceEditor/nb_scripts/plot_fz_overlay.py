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

# %% [markdown]
# # F vs Z overlay — viable zWLC fits (20260710)
#
# Overlays raw F(pN) vs Z(um) scatter, with the stored zWLC fit curve, for every
# ROI marked `selected` whose stored `R2 zWLC` exceeds a viability threshold.
#
# Selection and fit parameters (`L (bp)`, `A (nm)`, `S (pN)`, `Zbead (um)`, `R2 zWLC`)
# are read directly from each file's `/Parameters` table — nothing is recomputed.

# %%
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from trace_editor.core.h5_data import H5Data
from trace_editor.core.force_spectroscopy import zWLC, ExponentialDrift
from trace_editor.gui.main_window import _PALETTE_COLORS

STYLE_PATH = Path("../src/trace_editor/gui/traceeditor.mplstyle")
plt.style.use(STYLE_PATH)

FIT_COLOR = "black"
MEAN_COLOR = _PALETTE_COLORS["blue"]
BAND_COLOR = _PALETTE_COLORS["blue"]

Y_LIM = (0.4, 2.2)  # Z (um) axis limits, shared across all Z plots
AXLABEL_SIZE = 16
TICKLABEL_SIZE = 13
LEGEND_SIZE = 12

OUT_DIR = Path("../figures")
OUT_DIR.mkdir(exist_ok=True)

# %matplotlib inline

# %%
# --- configure here ---
date = 20260710 #lambda DNA
date = 20260629 #Giorgia's constructs
FOLDER = Path("/Users/barton/My Drive/data") / str(date)

R2_MIN = 0.9   # viability threshold on stored "R2 zWLC"
# ----------------------

FILES = sorted(FOLDER.glob("data_*.h5"))
print(f"{len(FILES)} h5 files in {FOLDER}")

# %%
# Discover viable traces: selected == 1 and R2 zWLC > R2_MIN, across all files.
viable = []  # list of dicts: path, roi, r2, L, A, S, zbead

for path in FILES:
    h5 = H5Data(path)
    for roi in h5.list_rois():
        if h5.get_selection(roi) != 1.0:
            continue
        r2 = h5.get_parameter("R2 full zWLC", roi)
        if not np.isfinite(r2) or r2 <= R2_MIN:
            continue
        viable.append(dict(
            path=path,
            roi=roi,
            r2=r2,
            L=h5.get_parameter("L (bp)", roi),
            A=h5.get_parameter("A (nm)", roi),
            S=h5.get_parameter("S (pN)", roi),
            zbead=h5.get_parameter("Zbead (um)", roi),
            drift_Z0=h5.get_parameter("Z0 (um)", roi, default=0.0),
            drift_dZ_dt=h5.get_parameter("dZ_dt (nm_s)", roi, default=0.0),
            drift_tau=h5.get_parameter("tau (s)", roi, default=np.inf),
        ))

print(f"{len(viable)} viable traces (R2 zWLC > {R2_MIN})")
for v in viable:
    print(f"  {v['path'].name}  roi={v['roi']}  R2={v['r2']:.3f}")

# %%
SMOOTH_WINDOW_FZ = 41  # samples, rolling-average window applied to Z before F-Z scatter

fig, ax = plt.subplots(figsize=(9, 7))

cmap = plt.get_cmap("tab20", max(len(viable), 1))

F_all = []
Z_traces = []  # (F, Z_smooth) per viable trace, for the average-trace plot later

for i, v in enumerate(viable):
    h5 = H5Data(v["path"])
    F = h5.get_shared("F (pN)")
    T = h5.get_shared("Time (s)")
    Z_raw = h5.get_trace("Z (um)", v["roi"])
    F_all.append(F)

    drift = ExponentialDrift(T, Z0=v["drift_Z0"], dZ_dt=v["drift_dZ_dt"], tau=v["drift_tau"])
    # Drift-corrected bead Z, then shifted into the WLC chain-extension frame
    # (Zbead maps the bead's drift-corrected zero to the chain's folded/near-zero state).
    Z = (Z_raw - drift) - v["zbead"]
    Z_smooth = pd.Series(Z).rolling(SMOOTH_WINDOW_FZ, center=True, min_periods=1).mean().to_numpy()
    Z_traces.append(Z_smooth)

    color = cmap(i)
    label = f"{v['path'].stem} roi{v['roi']} (R2={v['r2']:.2f})"

    ax.scatter(F, Z_smooth, s=6, alpha=0.45, color=color, edgecolors="none", label=label)

# zWLC model curves overlap across traces (same L/A/S) — draw once on top.
F_fit = np.linspace(max(np.concatenate(F_all).min(), 1e-3), np.concatenate(F_all).max(), 200)
Z_fit = zWLC(F_fit, L=5833, A=viable[0]["A"], S=viable[0]["S"], z_bead=0)  # L=viable[0]["L"]
ax.plot(F_fit, Z_fit, color=FIT_COLOR, linewidth=2.5, label="zWLC fit", zorder=100)

ax.set_xlabel("F (pN)", fontsize=AXLABEL_SIZE)
ax.set_ylabel("Z (µm)", fontsize=AXLABEL_SIZE)
ax.set_ylim(*Y_LIM)
ax.tick_params(labelsize=TICKLABEL_SIZE)
ax.legend(fontsize=LEGEND_SIZE, ncol=2, loc="lower right", frameon=True, framealpha=0.9)

fig.tight_layout()
fig.savefig(OUT_DIR / "fz_overlay_per_trace.pdf")
plt.show()

# %%
# Z vs T, same traces, drift- and Zbead-corrected (WLC reference frame).
SMOOTH_WINDOW = 40  # samples, rolling-average window for the Z(t) traces

fig, ax = plt.subplots(figsize=(9, 5))

cmap = plt.get_cmap("tab20", max(len(viable), 1))

T_traces = []  # per-trace (T, Z_smooth), for the average-trace plot later

for i, v in enumerate(viable):
    h5 = H5Data(v["path"])
    T = h5.get_shared("Time (s)")
    Z_raw = h5.get_trace("Z (um)", v["roi"])

    drift = ExponentialDrift(T, Z0=v["drift_Z0"], dZ_dt=v["drift_dZ_dt"], tau=v["drift_tau"])
    Z = (Z_raw - drift) - v["zbead"]
    Z_smooth = pd.Series(Z).rolling(SMOOTH_WINDOW, center=True, min_periods=1).mean().to_numpy()
    T_traces.append((T, Z_smooth))

    color = cmap(i)
    label = f"{v['path'].stem} roi{v['roi']} (R2={v['r2']:.2f})"

    ax.plot(T, Z_smooth, color=color, linewidth=1.0, alpha=0.85, label=label, zorder=i)

# zWLC-predicted Z(t) = zWLC(F(t)) — theory curve driven by the force protocol, on top.
h5_ref = H5Data(viable[0]["path"])
T_ref = h5_ref.get_shared("Time (s)")
F_ref = h5_ref.get_shared("F (pN)")
Z_wlc_t = zWLC(np.clip(F_ref, 1e-3, None), L=viable[0]["L"], A=viable[0]["A"], S=viable[0]["S"], z_bead=0)
ax.plot(T_ref, Z_wlc_t, color=FIT_COLOR, linewidth=2.5, label="zWLC fit: Z(F(t))", zorder=100)

ax.set_xlabel("Time (s)", fontsize=AXLABEL_SIZE)
ax.set_ylabel("Z (µm)", fontsize=AXLABEL_SIZE)
ax.set_xlim(0, T_ref[-1])
ax.set_ylim(*Y_LIM)
ax.tick_params(labelsize=TICKLABEL_SIZE)
# ax.legend(fontsize=LEGEND_SIZE, ncol=2, loc="best", frameon=True, framealpha=0.9)

fig.tight_layout()
fig.savefig(OUT_DIR / "z_vs_t_per_trace.pdf")
plt.show()

# %% [markdown]
# ## Average trace
#
# Traces run for different durations/sample counts, so each Z(t) is first
# resampled onto a common grid of *fractional* protocol time (0→1) before
# averaging. Shaded band = ±1 SD across traces at each point.

# %%
N_GRID = 300  # points in the common fractional-time grid
frac_grid = np.linspace(0, 1, N_GRID)

Z_resampled = np.empty((len(T_traces), N_GRID))
for i, (T, Z_smooth) in enumerate(T_traces):
    frac_t = (T - T[0]) / (T[-1] - T[0])
    Z_resampled[i] = np.interp(frac_grid, frac_t, Z_smooth)

Z_mean = Z_resampled.mean(axis=0)
Z_std = Z_resampled.std(axis=0)

# Use the longest-duration trace as the time axis for display.
i_ref = max(range(len(T_traces)), key=lambda i: T_traces[i][0][-1] - T_traces[i][0][0])
T_ref_avg = T_traces[i_ref][0]
T_display = T_ref_avg[0] + frac_grid * (T_ref_avg[-1] - T_ref_avg[0])

# zWLC-predicted Z(t) = zWLC(F(t)), resampled onto the same fractional grid.
F_ref_interp = np.interp(frac_grid, (T_ref - T_ref[0]) / (T_ref[-1] - T_ref[0]), F_ref)
Z_wlc_avg = zWLC(np.clip(F_ref_interp, 1e-3, None), L=viable[0]["L"], A=viable[0]["A"], S=viable[0]["S"], z_bead=0)

fig, ax = plt.subplots(figsize=(9, 5))

ax.fill_between(T_display, Z_mean - Z_std, Z_mean + Z_std,
                 color=BAND_COLOR, alpha=0.25, linewidth=0, label="±1 SD")
ax.plot(T_display, Z_mean, color=MEAN_COLOR, linewidth=2.5, label=f"Mean (n={len(T_traces)})")
ax.plot(T_display, Z_wlc_avg, color=FIT_COLOR, linewidth=2.5, label="zWLC fit: Z(F(t))", zorder=100)

ax.set_xlabel("Time (s)", fontsize=AXLABEL_SIZE)
ax.set_ylabel("Z (µm)", fontsize=AXLABEL_SIZE)
ax.set_xlim(0, T_display[-1])
ax.set_ylim(*Y_LIM)
ax.tick_params(labelsize=TICKLABEL_SIZE)
ax.legend(fontsize=LEGEND_SIZE, loc="best", frameon=True, framealpha=0.9)

fig.tight_layout()
fig.savefig(OUT_DIR / "z_vs_t_average.pdf")
plt.show()

# %% [markdown]
# ## Average trace — F vs Z
#
# Same averaged Z(t), now plotted against F (via the reference trace's F(t),
# resampled onto the same fractional grid) — mean ± 1 SD band, with the zWLC
# fit curve on top.

# %%
sort_idx = np.argsort(F_ref_interp)
F_sorted = F_ref_interp[sort_idx]
Z_mean_sorted = Z_mean[sort_idx]
Z_std_sorted = Z_std[sort_idx]

F_fit_avg = np.linspace(max(F_sorted.min(), 1e-3), F_sorted.max(), 200)
Z_fit_avg = zWLC(F_fit_avg, L=viable[0]["L"], A=viable[0]["A"], S=viable[0]["S"], z_bead=0)

fig, ax = plt.subplots(figsize=(9, 7))

ax.fill_between(F_sorted, Z_mean_sorted - Z_std_sorted, Z_mean_sorted + Z_std_sorted,
                 color=BAND_COLOR, alpha=0.25, linewidth=0, label="±1 SD")
ax.plot(F_sorted, Z_mean_sorted, color=MEAN_COLOR, linewidth=2.5, label=f"Mean (n={len(T_traces)})")
ax.plot(F_fit_avg, Z_fit_avg, color=FIT_COLOR, linewidth=2.5, label="zWLC fit", zorder=100)

ax.set_xlabel("F (pN)", fontsize=AXLABEL_SIZE)
ax.set_ylabel("Z (µm)", fontsize=AXLABEL_SIZE)
ax.set_ylim(*Y_LIM)
ax.tick_params(labelsize=TICKLABEL_SIZE)
ax.legend(fontsize=LEGEND_SIZE, loc="lower right", frameon=True, framealpha=0.9)

fig.tight_layout()
fig.savefig(OUT_DIR / "fz_average.pdf")
plt.show()
