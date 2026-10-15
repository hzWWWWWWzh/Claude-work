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
import sys

sys.path.insert(0, str(Path.cwd()))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trace_editor.core.h5_data import H5Data
from trace_editor.core.force_spectroscopy import zWLC, zFJC, ExponentialDrift
from trace_editor.gui.main_window import _PALETTE_COLORS

# %matplotlib inline

# %%
ROOT = Path("/Users/barton/1NAll/TraceEditor/hairpin/unzipping_analysis_package/")

HDF_DIR = ROOT / "data" / "hdf"          # raw TraceEditor .hdf files
FILES = ["data_015","data_016", "data_019", "data_020"]   # traces to analyse + merge
LABEL = "7"                 # bead label inside the .hdf

SEQUENCE_DNA = ROOT / "data" / "sequence" / "unzippingpartYT0309_49.dna"

SMOOTH_WINDOW_FZ = 50

# dsDNA handle contour length flanking the hairpin (NOT the 880 bp hairpin itself):
# 1.25 kb "Handle Right" + 1.05 kb "Handle Left" via the 3-way-junction adapter.
L_HANDLES_BP = 1250 + 1050

L_HAIRPIN_BP = 880  # fully unzipped hairpin length

# data_015 is ramped to ~20 pN, guaranteeing full unzip; used only to find the
# full-open anchor for rescaling raw N_unzip onto the L_HAIRPIN_BP scale.
ANCHOR_TRACE_NAME = "data_015"
OPEN_FORCE_MIN_PN = 13.0  # high-force cut used to find the full-open anchor


def baseline_offset_um(F, Z_meas_um, L_handles_bp, A, S, low_force_quantile=0.1):
    """
    Empirical Z offset (um) that anchors the closed-hairpin state (lowest-force
    points) to zero unzipped bp, since the drift-corrected Z=0 is an arbitrary
    bead position, not the physical zero-extension point of the WLC model.

    :param F: force (pN)
    :param Z_meas_um: measured, drift/Zbead-corrected extension (um)
    :param L_handles_bp: dsDNA handle contour length (bp)
    :param A: persistence length (nm)
    :param S: stretch modulus (pN)
    :param low_force_quantile: fraction of lowest-force points used to
        estimate the closed-hairpin baseline
    :return: offset (um) to add to Z_meas_um before calling unzipped_bp
    """
    threshold = np.nanquantile(F, low_force_quantile)
    baseline = F <= threshold
    z_handles_um = zWLC(F[baseline], L=L_handles_bp, A=A, S=S)
    return np.nanmean(z_handles_um - Z_meas_um[baseline])


def unzipped_bp(F, Z_meas_um, L_handles_bp, A, S):
    """
    Number of unzipped hairpin bp at each (F, Z) point.

    Z_meas is modeled as the extensible-WLC extension of the fixed dsDNA
    handles plus the FJC extension of the released single-stranded DNA
    (2 nt released per unzipped bp). Solving for the ssDNA contour length
    and dividing by the FJC extension-per-nt gives the unzipped bp count.

    :param F: force (pN)
    :param Z_meas_um: measured, drift/Zbead-corrected extension (um),
        already calibrated so the closed-hairpin state is at zero unzipped bp
        (see baseline_offset_um)
    :param L_handles_bp: dsDNA handle contour length (bp)
    :param A: persistence length (nm)
    :param S: stretch modulus (pN)
    :return: unzipped hairpin length (bp)
    """
    z_handles_um = zWLC(F, L=L_handles_bp, A=A, S=S)
    z_ss_nm = (Z_meas_um - z_handles_um) * 1000
    z_per_nt_nm = zFJC(F, L=1)
    N_unzip_bp = z_ss_nm / (2 * z_per_nt_nm)
    return N_unzip_bp


# %%
viable = []
for path in FILES:
    full_path = HDF_DIR/f"{path}.h5"
    h5 = H5Data(full_path)
    for roi in h5.list_rois():
        # if h5.get_selection(roi) != 1.0:
        #         continue
        if roi == '007':
            viable.append(dict(path=full_path, name=path, roi=roi, h5=h5))
viable

# %%
def scatter_with_hist(fig, gs, row, col, raw_by_trace, smooth_by_trace, ylabel):
    """
    Plot concatenated per-trace raw (grey) + smoothed (colored) series, with
    a shared-y histogram of the smoothed values.
    """
    trace_ax = fig.add_subplot(gs[row, col])
    hist_ax = fig.add_subplot(gs[row, col + 1], sharey=trace_ax)

    start = 0
    for (name, raw), (_, smooth) in zip(raw_by_trace, smooth_by_trace):
        x = np.arange(start, start + len(raw))
        trace_ax.plot(x, raw, '.', color='grey', alpha=0.3, markersize=2)
        line, = trace_ax.plot(x, smooth, '-', label=name)
        hist_ax.hist(smooth, bins=40, orientation="horizontal", color=line.get_color(), alpha=0.35)
        start += len(raw)

    trace_ax.set(xlabel="Concatenated sample", ylabel=ylabel)
    trace_ax.grid(alpha=0.3)
    trace_ax.set_xlim(0, start)
    hist_ax.set(xlabel="Count")
    plt.setp(hist_ax.get_yticklabels(), visible=False)
    return trace_ax, hist_ax


# %%
# Z (raw and smoothed) and unzipped length (raw and smoothed), pooled across
# traces; smoothed N_unzip histogram reveals stalling positions as peaks.
Z_by_trace, Z_smooth_by_trace = [], []
N_unzip_by_trace, N_unzip_smooth_by_trace = [], []
F_by_trace = []

for trace in viable:
    h5, roi = trace["h5"], trace["roi"]
    F = h5.get_shared("F (pN)")
    T = h5.get_shared("Time (s)")
    Z_raw = h5.get_trace("Z (um)", roi)

    A = h5.get_parameter("A (nm)", roi, default=50.0)
    S = h5.get_parameter("S (pN)", roi, default=1000.0)
    zbead = h5.get_parameter("Zbead (um)", roi, default=0.0)
    drift = ExponentialDrift(
        T,
        Z0=h5.get_parameter("Z0 (um)", roi, default=0.0),
        dZ_dt=h5.get_parameter("dZ_dt (nm_s)", roi, default=0.0),
        tau=h5.get_parameter("tau (s)", roi, default=np.inf),
    )
    Z = (Z_raw - drift) - zbead
    Z += baseline_offset_um(F, Z, L_HANDLES_BP, A, S)
    Z_smooth = pd.Series(Z).rolling(SMOOTH_WINDOW_FZ, center=True, min_periods=1).mean().to_numpy()

    N_unzip = unzipped_bp(F, Z, L_HANDLES_BP, A, S)
    N_unzip_smooth = unzipped_bp(F, Z_smooth, L_HANDLES_BP, A, S)

    Z_by_trace.append((trace["name"], Z))
    Z_smooth_by_trace.append((trace["name"], Z_smooth))
    N_unzip_by_trace.append((trace["name"], N_unzip))
    N_unzip_smooth_by_trace.append((trace["name"], N_unzip_smooth))
    F_by_trace.append((trace["name"], F))

# Rescale N_unzip so the full-open anchor (from the guaranteed-full-unzip
# ANCHOR_TRACE_NAME) matches L_HAIRPIN_BP; the dWLC/FJC model parameters are
# only approximate, so N_unzip is proportional to bp but not exactly on scale.
# The anchor uses the heavily smoothed trace's high-force plateau, since
# smoothing already suppresses tracking noise there.
anchor_F, anchor_N_unzip_smooth = next(
    (F, N_unzip)
    for (name, F), (_, N_unzip) in zip(F_by_trace, N_unzip_smooth_by_trace)
    if name == ANCHOR_TRACE_NAME
)
open_raw = np.nanmax(anchor_N_unzip_smooth[anchor_F >= OPEN_FORCE_MIN_PN])
scale = L_HAIRPIN_BP / open_raw
print(f"full-open anchor: {open_raw:.1f} bp (raw) -> scale factor {scale:.4f}")

N_unzip_by_trace = [(name, N_unzip * scale) for name, N_unzip in N_unzip_by_trace]
N_unzip_smooth_by_trace = [(name, N_unzip * scale) for name, N_unzip in N_unzip_smooth_by_trace]

fig = plt.figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], wspace=0.05, hspace=0.3)

top_ax, _ = scatter_with_hist(fig, gs, 0, 0, Z_by_trace, Z_smooth_by_trace, "Z (µm)")

force_ax = top_ax.twinx()
start = 0
for name, F in F_by_trace:
    x = np.arange(start, start + len(F))
    force_ax.plot(x, F, '-', color='grey', alpha=0.25, linewidth=1, zorder=0)
    start += len(F)
force_ax.set(ylabel="F (pN)")
force_ax.set_zorder(top_ax.get_zorder() - 1)
top_ax.patch.set_visible(False)
bottom_ax, bottom_hist_ax = scatter_with_hist(
    fig, gs, 1, 0, N_unzip_by_trace, N_unzip_smooth_by_trace, "Unzipped length (bp)"
)
bottom_ax.axhline(0, color='k', linestyle='--', linewidth=1, label="fully zipped")
bottom_ax.axhline(L_HAIRPIN_BP, color='r', linestyle='--', linewidth=1, label="fully unzipped")
bottom_hist_ax.axhline(0, color='k', linestyle='--', linewidth=1)
bottom_hist_ax.axhline(L_HAIRPIN_BP, color='r', linestyle='--', linewidth=1)
# top_ax.legend()
# bottom_ax.legend()
plt.show()

# %%
# TEMP: overlay Z_smooth and N_unzip_smooth (twin axes) to sanity-check that
# the WLC/FJC unzipping conversion is doing more than a rescale of Z.
fig, axes = plt.subplots(len(viable), 1, figsize=(10, 3 * len(viable)))
axes = np.atleast_1d(axes)

for ax, trace, (_, Z_smooth), (_, N_unzip_smooth) in zip(
    axes, viable, Z_smooth_by_trace, N_unzip_smooth_by_trace
):
    T = trace["h5"].get_shared("Time (s)")
    l1, = ax.plot(T, Z_smooth, color="C0", label="Z smoothed (µm)")
    ax.set(title=trace["name"], xlabel="Time (s)", ylabel="Z (µm)")
    ax.yaxis.label.set_color(l1.get_color())

    ax2 = ax.twinx()
    l2, = ax2.plot(T, N_unzip_smooth, color="C1", label="N_unzip smoothed (bp)")
    ax2.set(ylabel="Unzipped length (bp)")
    ax2.yaxis.label.set_color(l2.get_color())

    ax.legend(handles=[l1, l2], loc="upper left")
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# %%
# Simulate expected unzipping behavior from the hairpin sequence, since the
# force channel calibration is uncertain and this gives an independent
# ground truth to compare measured traces against.
import struct


def read_snapgene_sequence(path):
    """
    Extract the raw base sequence from a SnapGene .dna file (block type 0:
    1 topology byte + ASCII bases). SnapGene stores the file as a series of
    (1-byte type, 4-byte big-endian length, payload) blocks.

    :param path: path to a .dna file
    :return: sequence (uppercased str)
    """
    data = Path(path).read_bytes()
    pos = 0
    while pos + 5 <= len(data):
        block_type = data[pos]
        length = struct.unpack(">I", data[pos + 1:pos + 5])[0]
        payload_start = pos + 5
        if block_type == 0:
            return data[payload_start + 1:payload_start + length].decode("ascii").upper()
        pos = payload_start + length
    raise ValueError(f"{path}: no sequence block (type 0) found")


HAIRPIN_SEQUENCE = read_snapgene_sequence(SEQUENCE_DNA)[::1]
print(f"hairpin sequence: {len(HAIRPIN_SEQUENCE)} bp")

# %%
# SantaLucia (1998) unified nearest-neighbor DeltaG_37 parameters (kcal/mol,
# 1M NaCl); literature values, standard in the unzipping-DNA field, not
# specific to any codebase.
NN_DELTA_G_KCAL_MOL = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
KCAL_PER_MOL_PER_KT = 0.593  # at ~298 K

SS_CONTOUR_NM_PER_BASE = 0.58
SS_KUHN_NM = 1.55
SS_STRETCH_PN = 800.0


def nn_opening_cost_kt(sequence):
    """
    Free energy cost (kT) to open each successive base-pair step, from
    nearest-neighbor duplex stability. NN_DELTA_G is the (negative) free
    energy of duplex FORMATION for each dinucleotide step, so the cost to
    break it is the negated, kT-converted value.

    :param sequence: dsDNA sequence (5'->3', uppercase ACGT), read from the
        end that unzips first towards the end that unzips last
    :return: array of length len(sequence)-1, opening cost (kT) per step
    """
    costs = []
    for i in range(len(sequence) - 1):
        step = sequence[i:i + 2]
        dg_kcal_mol = NN_DELTA_G_KCAL_MOL.get(step)
        if dg_kcal_mol is None:
            dg_kcal_mol = np.nan
        costs.append(-dg_kcal_mol / KCAL_PER_MOL_PER_KT)
    return np.array(costs)


def _langevin(x):
    x = np.atleast_1d(x).astype(float)
    out = np.empty_like(x)
    small = np.abs(x) < 1e-4
    out[small] = x[small] / 3.0
    xs = x[~small]
    out[~small] = np.cosh(xs) / np.sinh(xs) - 1.0 / xs
    return out.item() if out.size == 1 else out


def mechanical_work_per_step_kt(force_pn):
    """
    Mechanical work gained (kT) by releasing one bp of the hairpin stem into
    2 nt of ssDNA, stretched at the given force. The hairpin stem is not
    under tension while zipped (it is not part of the force-bearing path
    between the two dsDNA handles), so opening a bp adds pure ssDNA
    extension with no dsDNA term to subtract.

    :param force_pn: force (pN)
    :return: work per opened bp (kT)
    """
    f = max(force_pn, 1e-6)
    x = f * SS_KUHN_NM / 4.114
    ss_ext_nm = SS_CONTOUR_NM_PER_BASE * (_langevin(x) + f / SS_STRETCH_PN)
    delta_ext_nm = 2 * ss_ext_nm
    return f * delta_ext_nm / 4.114


def free_energy_landscape(sequence, force_pn):
    """
    Free energy G(N) (kT, relative to the fully closed state N=0) as a
    function of unzipped length N, at a single force: sequence-dependent NN
    opening cost minus mechanical work gained from being pulled at that
    force. G(0) = 0 always -- the closed state is a fixed physical
    reference, not the landscape's own (force-dependent) global minimum,
    so this stays comparable across different forces.

    :param sequence: dsDNA sequence of the hairpin (5'->3', unzips in order)
    :param force_pn: force (pN)
    :return: G(N) array, length len(sequence), with G(0) = 0
    """
    cum_cost_kt = np.concatenate([[0.0], np.cumsum(nn_opening_cost_kt(sequence))])
    N = np.arange(len(cum_cost_kt))
    return cum_cost_kt - N * mechanical_work_per_step_kt(force_pn)


def simulate_unzip(sequence, forces_pn):
    """
    Equilibrium unzipped-length distribution at each force, from a
    partition-function model over unzipped length N: opening the first N bp
    costs sum(opening_cost[:N]) in sequence-dependent NN free energy, and
    gains N * mechanical_work_per_step(F) from being pulled at force F.

    :param sequence: dsDNA sequence of the hairpin (5'->3', unzips in order)
    :param forces_pn: forces (pN) to evaluate
    :return: (mean_N_unzip, P) where P is (len(forces_pn), len(sequence))
        the equilibrium occupancy of each unzipped-length state at each
        force (N ranges 0..len(sequence)-1)
    """
    mean_N = np.empty(len(forces_pn))
    P = np.empty((len(forces_pn), len(sequence)))
    for i, f in enumerate(forces_pn):
        G = free_energy_landscape(sequence, f)
        N = np.arange(len(G))
        p = np.exp(-G)
        p /= p.sum()
        P[i] = p
        mean_N[i] = np.sum(N * p)
    return mean_N, P


def kramers_dwell_time_s(sequence, force_pn, k0_per_s=1e6, smooth_bp=25):
    """
    Kramers escape-rate dwell time (s) at every local free-energy minimum of
    G(N) at a given force. The barrier confining a minimum is the height up
    to the nearest local maximum on either side (or the landscape edge, if
    no local maximum exists on that side) -- G(N) has no separate
    transition-state energy, so this treats the 1D landscape's own
    roughness as the barrier, standard for a Kramers estimate on a rugged
    free-energy landscape.

    G(N) is coarse-grained over smooth_bp before finding minima/maxima: at
    single-bp resolution the landscape is dominated by sequence ruggedness
    (real, but finer than a real pause -- every AT/GC step creates a
    trivial ~1 kT wiggle), which floods the result with spurious minima.

    :param sequence: dsDNA sequence of the hairpin (5'->3', unzips in order)
    :param force_pn: force (pN)
    :param k0_per_s: attempt frequency (1/s); sets the absolute timescale
        only -- relative dwell times between minima do not depend on it
    :param smooth_bp: coarse-graining window (bp) applied to G(N) before
        finding minima/maxima
    :return: (N_minima, dwell_time_s) arrays for each local minimum found
    """
    G_raw = free_energy_landscape(sequence, force_pn)
    G = pd.Series(G_raw).rolling(smooth_bp, center=True, min_periods=1).mean().to_numpy()
    n = len(G)

    is_min = np.zeros(n, dtype=bool)
    is_min[1:-1] = (G[1:-1] <= G[:-2]) & (G[1:-1] <= G[2:])
    minima = np.where(is_min)[0]

    is_max = np.zeros(n, dtype=bool)
    is_max[1:-1] = (G[1:-1] >= G[:-2]) & (G[1:-1] >= G[2:])
    maxima = np.where(is_max)[0]

    dwell_s = np.empty(len(minima))
    for j, N in enumerate(minima):
        left_maxima = maxima[maxima < N]
        right_maxima = maxima[maxima > N]
        barrier_left = G[left_maxima[-1]] if left_maxima.size else G[0]
        barrier_right = G[right_maxima[0]] if right_maxima.size else G[-1]
        rate = k0_per_s * (np.exp(-(barrier_left - G[N])) + np.exp(-(barrier_right - G[N])))
        dwell_s[j] = 1.0 / rate
    return minima, dwell_s

# %%
# Simulated <N_unzip>(F) curve and equilibrium occupancy map, plus a
# force-ramp-order dwell-time proxy: forces where P(N,F) has multiple
# comparably-populated N are where the real trace should stall/hop.
sim_forces = np.linspace(0.5, 20, 300)
sim_mean_N, sim_P = simulate_unzip(HAIRPIN_SEQUENCE[::1], sim_forces)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(sim_forces, sim_mean_N, color="k")
axes[0].set(xlabel="F (pN)", ylabel="Simulated <unzipped length> (bp)")
axes[0].grid(alpha=0.3)

im = axes[1].imshow(
    sim_P.T, origin="lower", aspect="auto",
    extent=[sim_forces.min(), sim_forces.max(), 0, sim_P.shape[1]],
    cmap="viridis",
)
axes[1].set(xlabel="F (pN)", ylabel="Unzipped length (bp)", title="Equilibrium occupancy P(N | F)")

plt.tight_layout()
plt.show()

# %%
# Kramers dwell time at each local free-energy minimum, scanned across the
# transition region where local minima actually exist (outside it -- e.g.
# fully closed at low F or fully open at high F -- the landscape is
# monotonic and there is nothing to stall at). k0 = 1e6 /s is a typical
# literature order-of-magnitude estimate for single-bp opening attempt
# rates -- only the RELATIVE dwell times between minima are meaningful,
# not the absolute seconds.
DWELL_SCAN_FORCES_PN = np.arange(10.0, 11.5, 0.2)

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
cmap = plt.get_cmap("viridis", len(DWELL_SCAN_FORCES_PN))

for i, f in enumerate(DWELL_SCAN_FORCES_PN):
    G_smooth = pd.Series(free_energy_landscape(HAIRPIN_SEQUENCE, f)).rolling(
        25, center=True, min_periods=1
    ).mean().to_numpy()
    minima_N, dwell_s = kramers_dwell_time_s(HAIRPIN_SEQUENCE, f, k0_per_s=1e6, smooth_bp=25)
    color = cmap(i)
    axes[0].plot(np.arange(len(G_smooth)), G_smooth, color=color, linewidth=1, label=f"{f:.1f} pN")
    axes[0].plot(minima_N, G_smooth[minima_N], 'o', color=color)
    axes[1].plot(minima_N, dwell_s, 'o-', color=color, label=f"{f:.1f} pN")

axes[0].set(ylabel="G(N) smoothed (kT)", title="Coarse-grained free energy landscape across the transition region")
axes[0].legend(fontsize=8, ncol=3)
axes[0].grid(alpha=0.3)

axes[1].set(xlabel="Unzipped length N (bp)", ylabel="Dwell time (s)", yscale="log")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

# %%
# Overlay simulated <N_unzip>(F) on top of measured N_unzip_smooth(F) per
# trace, to check whether the sequence-predicted dwell/stall forces line up
# with the measured data independent of the force-channel calibration issue.
fig, axes = plt.subplots(len(viable), 1, figsize=(10, 3.5 * len(viable)))
axes = np.atleast_1d(axes)

for ax, trace, (_, F), (_, N_unzip_smooth) in zip(axes, viable, F_by_trace, N_unzip_smooth_by_trace):
    ax.plot(F, N_unzip_smooth, '.', markersize=2, alpha=0.4, color="C0", label="measured (smoothed)")
    ax.plot(sim_forces, sim_mean_N, color="C3", label="simulated <N_unzip>(F)")
    ax.set(title=trace["name"], xlabel="F (pN)", ylabel="Unzipped length (bp)")
    ax.legend()
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()
