import warnings

warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=RuntimeWarning)

import TraceIO as tio
from pathlib import Path
from CustomIcecream import ic
from matplotlib import pyplot as plt
import ForceSpectroscopy as fs
import numpy as np
import pandas as pd
from scipy.optimize import approx_fprime, minimize
from scipy.ndimage import median_filter
from scipy.signal import find_peaks


def get_zbead(filename, label, save=False, show=False):
    ic(filename)
    data = tio.hdf_data(filename, label)
    # print(data.list_channels())
    # ic(data.list_parameters())

    z = data.traces["Z (um)"]
    f = data.traces["F (pN)"]
    t = data.traces["Time (s)"]
    s = data.traces["Shift (mm)"]
    l = data.get_parameter("L (bp)")

    selection = s < np.percentile(s, 20)  # discard very low forces

    if not selection.any():
        print("Warning: No data points selected for analysis!")
        selection = np.ones_like(s, dtype=bool)  # Use all points as fallback

    z -= data.calculate_drift()
    z_wlc = fs.zWLC(f, l)

    k = np.diff(f) / np.diff(z_wlc)
    k = np.concatenate(([k[0]], k))
    # replace nan by nearest value
    k = pd.Series(k).interpolate(method="nearest").bfill().ffill().values

    # replace zero or negative k by nearest positive value
    mask = k <= 0
    if mask.any():
        k_positive = k.copy()
        k_positive[mask] = np.nan
        k[mask] = (
            pd.Series(k_positive)
            .interpolate(method="nearest")
            .bfill()
            .ffill()
            .values[mask]
        )

    kT = 4.1e-3  # pN um at room temperature
    sigma_z = np.sqrt(kT / k)

    tracking_error = 10e-3  # in um
    sigma_z = np.sqrt(sigma_z**2 + tracking_error**2)

    def neg_goodness(z_offset):
        z_scores = (z[selection] - (z_wlc[selection] + z_offset)) / sigma_z[selection]
        return -np.exp(
            -0.5 * z_scores**2
        ).mean()  # negative because we want to maximize

    # Initial guess
    z_guess = np.percentile(z, 90) - np.percentile(z_wlc, 90)

    result = minimize(
        neg_goodness,
        x0=z_guess,
        bounds=[(z_guess - 3, z_guess + 3)],
        method="L-BFGS-B",
        options={"ftol": 1e-9},
    )  # increase precision
    z_bead_opt = result.x[0]

    # Estimate error using Hessian (curvature at the minimum)
    eps = np.sqrt(np.finfo(float).eps)
    hess_approx = approx_fprime(
        [z_bead_opt], lambda x: approx_fprime(x, neg_goodness, eps), eps
    )[0]
    z_bead_error = np.sqrt(1.0 / abs(hess_approx)) if hess_approx != 0 else np.nan

    z_bead_string = tio.round_significance(z_bead_opt, z_bead_error)

    z_scores = (z - (z_wlc + z_bead_opt)) / sigma_z
    goodness = np.exp(-0.5 * z_scores**2)

    if save:
        data.set_parameter("Zbead (um)", z_bead_opt, error=z_bead_error, type="fit")
        data.save()

    if show:
        # plot z vs time and z_scores vs time in two plots below each other
        plt.figure(figsize=(10, 8))
        plt.subplot(2, 1, 1)
        plt.plot(t, goodness, ".", color="r", markersize=5, alpha=0.1)
        plt.xlabel("Time (s)")
        plt.ylabel("Goodness")
        plt.subplot(2, 1, 2)
        plt.fill_between(t, z_wlc - sigma_z, z_wlc + sigma_z, color="gray", alpha=0.3)
        plt.plot(t, z - z_bead_opt, ".", color="gray", markersize=5, alpha=0.1)
        plt.plot(
            t[selection],
            z[selection] - z_bead_opt,
            ".",
            color="b",
            markersize=5,
            alpha=0.1,
        )
        plt.xlabel("Time (s)")
        plt.ylabel("Z (um)")
        plt.plot(t, z_wlc, color="k")
        plt.title(f"z_bead = {z_bead_string} (um)")
        plt.ylim(bottom=0, top=np.max(z_wlc) * 1.1)
        plt.hlines(-z_bead_opt, color="k", linestyle="--", xmin=t.min(), xmax=t.max())
        plt.tight_layout()
        plt.show()
    return z_bead_opt, z_bead_error


def get_steps(filename, label, n=30, l2=80, show=False):
    data = tio.hdf_data(filename, label)
    z = data.traces["Z (um)"] - data.calculate_drift()
    f = data.traces["F (pN)"]
    t = data.traces["Time (s)"]
    l = data.get_parameter("L (bp)")

    z -= get_zbead(filename, label, save=True, show=False)[0]
    z_wlc = fs.zWLC(f, l)

    extensions = [fs.zWLC(f, l - i * l2) for i in range(n + 1)]
    diffs = np.asarray([abs(z - e) for e in extensions])

    best_l = np.argmin(diffs, axis=0)
    best_l = median_filter(best_l, size=50)  # take median to get a single best fit

    z_fit = np.array([extensions[l][i] for i, l in enumerate(best_l)])

    if show:
        plt.plot(z, f, ".", color="gray", markersize=5, alpha=1)
        plt.plot(z_fit, f, "-", color="r", markersize=5, alpha=1)
        plt.xlabel("Z (um)")
        plt.ylabel("F (pN)")
        plt.show()

    return z_fit, best_l


def zFiber(
    x,
    L=6000,
    A=50,
    S=1000,
    n0=0,
    n4=0,
    n=0,
    k=0.1,
    NLD=1.2,
    dg_unwrap=12,
    dg1=17,
    dg2=5,
    dg3=50,
    NRL=197,
    L0=147,
    L1=89,
    L2=80,
    d=1,
    z_bead=0,
):
    """
    Compute force-dependent extension of chromatin fiber according to Meng2019NAR

    :param x: F (pN) # Force
    :param L: L (bp) #contour length (bp)
    :param A: A (nm) # persistence length
    :param S: S (pN) # Stretch modulus
    :param n0: n0 # number of nucleosomes that don not stack into a fiber
    :param n4: n4 # number of tetrasomes
    :param n: n # number of nucleosomes (incl. tetrasomes)
    :param k: k (pN_nm) # stiffness (pN/nm)
    :param NLD: NLD (nm) #Nucleosome Line Density (nm)
    :param dg_unwrap: Gwrap (kT) # wrapping energy for a free (non-stacked) nucleosome (kT)
    :param dg1: G1 (kT) #s tacking energy (kT)
    :param dg2: G2 (kT) # partial unwrapping energy (kT) of (l2-L1) basepairs
    :param dg3: G3 (kT) # full unwrapping energy (kT) of last L2 basepairs
    :param NRL: NRL (bp) # Nucleosome Repeat Length (bp)
    :param L0: L0 (bp) # wrapped amount of DNA (bp) in a wrapped, non-stacked nucleosome
    :param L1: L1 (bp) # wrapped amount of DNA (bp) in unstacked conformation
    :param L2: L2 (bp) # wrapped amount of DNA (bp) in singly wrapped conformation
    :param d: degeneracy # degeneracy (0 = cooperative, 1 = degenerate)
    :param z_bead: Zbead (um) # Height offset
    :return: Z (um) # fiber extension
    """

    x = np.asarray(x)

    # Compute extension and energy per stacked nucleosome and per basepair
    z0 = fs.zFJC(x, k=k, z0=NLD, L=25)
    g0 = fs.gHooke(x, k, NLD)

    # Compute the extension and energy of free DNA
    z_bp = fs.zWLC(x, 1, S=S, A=A) * 1000
    g_bp = fs.gWLC(x, 1, S=S, A=A) / fs.kT

    # Combine extension and energy of each states for all forces
    # thorder of states: 0 = stacked, 1 = partially unwrapped,
    # 2 = single wrapp, 3 = fully unwrapped
    z_state = np.array(
        [
            z0,
            z_bp * (NRL - L1),
            z_bp * (NRL - L2),
            z_bp * NRL,
        ]
    )
    g_state = np.array(
        [
            g0,
            g_bp * (NRL - L1) + dg1,
            g_bp * (NRL - L2) + dg1 + dg2,
            g_bp * NRL + dg1 + dg2 + dg3,
        ]
    )

    # subtract work done by force
    g_state = g_state - x[np.newaxis, :] * z_state / fs.kT
    g_min = np.min(g_state, axis=0)  # Avoid numerical overflow
    g_state -= g_min[np.newaxis, :]

    # compute Boltzmann weights
    P_state = np.exp(-g_state)
    sum = np.sum(P_state, axis=0)
    P_state /= sum[np.newaxis, :]

    if d < 1:  # adjust for cooperativity of first transition
        g_state = (1 - d) * n * g_state[:2]
        g_min = np.min(g_state, axis=0)  # Avoid numerical overflow
        g_state -= g_min[np.newaxis, :]
        P_state_coop = np.exp(-g_state)
        sum = np.sum(P_state_coop, axis=0)
        P_state_coop /= sum[np.newaxis, :]
        P_state[0] = P_state_coop[0]
        P_state[1] *= 0
        P_state[1] = 1 - np.sum(P_state, axis=0)

    # compute average extension
    z_fiber = np.sum(P_state * z_state, axis=0)
    z_fiber *= n - n4

    if n4 > 0:  # add tetrasomes
        p4 = np.exp(-g_state[2:])
        sum = np.sum(p4, axis=0)
        p4 /= sum[np.newaxis, :]
        z_fiber += np.sum(p4 * z_state[2:], axis=0) * n4

    z_fiber += z_bp * (L - n * NRL)  # add DNA handles

    return z_fiber * 0.001 + z_bead


def z_nucleosome(
    x,
    L=1000,
    A=50,
    S=1000,
    NRL=197,
    dg_unwrap=12,
    dg2=5,
    dg3=50,
    L0=147,
    L1=89,
    L2=80,
):
    """
    Compute force-dependent extension of chromatin fiber according to Meng2019NAR

    :param x: F (pN) # Force
    :param L: L (bp) #contour length (bp)
    :param A: A (nm) # persistence length
    :param S: S (pN) # Stretch modulus
    :param n0: n0 # number of nucleosomes that don not stack into a fiber
    :param n4: n4 # number of tetrasomes
    :param n: n # number of nucleosomes (incl. tetrasomes)
    :param k: k (pN_nm) # stiffness (pN/nm)
    :param NLD: NLD (nm) #Nucleosome Line Density (nm)
    :param dg_unstacked: Gwrap (kT) # wrapping energy for a free (non-stacked) nucleosome (kT)
    :param dg1: G1 (kT) #s tacking energy (kT)
    :param dg2: G2 (kT) # partial unwrapping energy (kT) of (l2-L1) basepairs
    :param dg3: G3 (kT) # full unwrapping energy (kT) of last L2 basepairs
    :param NRL: NRL (bp) # Nucleosome Repeat Length (bp)
    :param L0: L0 (bp) # wrapped amount of DNA (bp) in a wrapped, non-stacked nucleosome
    :param L1: L1 (bp) # wrapped amount of DNA (bp) in unstacked conformation
    :param L2: L2 (bp) # wrapped amount of DNA (bp) in singly wrapped conformation
    :param d: degeneracy # degeneracy (0 = cooperative, 1 = degenerate)
    :param z_bead: Zbead (um) # Height offset
    :return: Z (um) # fiber extension
    """

    x = np.asarray(x)
    # Compute the extension and energy of a single, unstacked nucleoseome
    # kinked DNA path implies shorter persistence length
    z_wrapped = fs.zWLC(x, NRL - 147, S=S, A=30)
    g_wrapped = fs.gWLC(x, NRL - 147, S=S, A=30) / fs.kT

    # Compute the extension and energy of free DNA
    z_bp = fs.zWLC(x, 1, S=S, A=A) * 1000
    g_bp = fs.gWLC(x, 1, S=S, A=A) / fs.kT

    # Combine extension and energy of each states for all forces
    # the order of states: 0 = fully wrapped, 1 = partially unwrapped,
    # 2 = single wrapp, 3 = fully unwrapped
    z_state = np.array(
        [
            z_wrapped,
            z_bp * (NRL - L1),
            z_bp * (NRL - L2),
            z_bp * NRL,
        ]
    )
    g_state = np.array(
        [
            g_wrapped,
            g_bp * (NRL - L1) + dg_unwrap,
            g_bp * (NRL - L2) + dg_unwrap + dg2,
            g_bp * NRL + dg_unwrap + dg2 + dg3,
        ]
    )

    # subtract work done by force
    g_state = g_state - x[np.newaxis, :] * z_state / fs.kT
    g_min = np.min(g_state, axis=0)  # Avoid numerical overflow
    g_state -= g_min[np.newaxis, :]

    # compute Boltzmann weights
    P_state = np.exp(-g_state)
    sum = np.sum(P_state, axis=0)
    P_state /= sum[np.newaxis, :]

    # compute average extension
    z_nucleosome = np.sum(P_state * z_state, axis=0)

    z_nucleosome += z_bp * (L - NRL)  # add DNA handles

    return z_nucleosome * 0.001


if __name__ == "__main__":
    label = "16"
    # Use absolute path relative to this script's location
    script_dir = Path(__file__).parent
    filename = script_dir / "data" / "data_008.hdf"

    if not filename.exists():
        print(f"Error: File {filename} not found!")
        print(f"Current directory: {Path.cwd()}")
        print(f"Script directory: {script_dir}")
        exit(1)

    # get_steps(filename, label, show=True)

    L = 1497

    force = np.logspace(-2, 1.8, 100)
    z = z_nucleosome(
        force,
        L=L,
        A=50,
        S=1000,
        NRL=197,
        dg_unwrap=19,
        dg2=5,
        dg3=50,
        L0=147,
        L1=89,
        L2=80,
    )
    plt.plot(z, force, "-", color="gray", markersize=5, alpha=1)
    plt.xlabel("Z (um)")
    plt.ylabel("F (pN)")
    # plt.yscale("log")
    plt.xlim(0, 1.4*L/3000)
    plt.show()
    

