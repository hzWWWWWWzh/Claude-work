import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=RuntimeWarning)

import TraceIO as tio
from pathlib import Path
from icecream import ic
from matplotlib import pyplot as plt
import ForceSpectroscopy as fs
import numpy as np
import pandas as pd
from scipy.optimize import approx_fprime, minimize
from scipy.ndimage import median_filter



def get_zbead(filename, label, save = False,show = False):
    ic(filename)
    data = tio.hdf_data(filename, label)
    # print(data.list_channels())
    # ic(data.list_parameters())

    z = data.traces['Z (um)']
    f = data.traces['F (pN)']
    t = data.traces['Time (s)']
    s = data.traces['Shift (mm)']
    l = data.get_parameter('L (bp)')

    selection = s < np.percentile(s, 20) # discard very low forces

    z -= data.calculate_drift()
    z_wlc = fs.zWLC(f, l)


    k = np.diff(f) / np.diff(z_wlc)
    k = np.concatenate(([k[0]], k))
    # replace nan by nearest value
    k = pd.Series(k).interpolate(method='nearest').bfill().ffill().values

    # replace zero or negative k by nearest positive value
    mask = k <= 0
    if mask.any():
        k_positive = k.copy()
        k_positive[mask] = np.nan
        k[mask] = pd.Series(k_positive).interpolate(method='nearest').bfill().ffill().values[mask]

    kT = 4.1e-3  # pN um at room temperature
    sigma_z = np.sqrt(kT/k)

    tracking_error = 10e-3  # in um
    sigma_z = np.sqrt(sigma_z**2 + tracking_error**2)



    def neg_goodness(z_offset):
        z_scores = (z[selection] - (z_wlc[selection] + z_offset)) / sigma_z[selection]
        return -np.exp(-0.5 * z_scores**2).mean()  # negative because we want to maximize

    # Initial guess
    z_guess = np.percentile(z, 90) - np.percentile(z_wlc, 90)

    result = minimize(neg_goodness, x0=z_guess, bounds=[(z_guess-3, z_guess+3)], method='L-BFGS-B',
                    options={'ftol': 1e-9})  # increase precision
    z_bead_opt = result.x[0]

    # Estimate error using Hessian (curvature at the minimum)
    eps = np.sqrt(np.finfo(float).eps)
    hess_approx = approx_fprime([z_bead_opt], lambda x: approx_fprime(x, neg_goodness, eps), eps)[0]
    z_bead_error = np.sqrt(1.0 / abs(hess_approx)) if hess_approx != 0 else np.nan

    z_bead_string = (tio.round_significance(z_bead_opt, z_bead_error))

    z_scores = (z - (z_wlc +z_bead_opt)) / sigma_z
    goodness = np.exp(-0.5 * z_scores**2)

    if save:
        data.set_parameter('Zbead (um)', z_bead_opt, error=z_bead_error, type='fit')
        data.save()

    if show:
        # plot z vs time and z_scores vs time in two plots below each other
        plt.figure(figsize=(10, 8))
        plt.subplot(2, 1, 1)
        plt.plot(t, goodness, '.', color='r', markersize=5, alpha=0.1)
        plt.xlabel('Time (s)')
        plt.ylabel('Goodness')
        plt.subplot(2, 1, 2)
        plt.fill_between(t, z_wlc - sigma_z, z_wlc + sigma_z, color='gray', alpha=0.3)
        plt.plot(t, z - z_bead_opt, '.', color='gray', markersize=5, alpha=0.1)
        plt.plot(t[selection], z[selection] - z_bead_opt, '.', color='b', markersize=5, alpha=0.1)
        plt.xlabel('Time (s)')
        plt.ylabel('Z (um)')
        plt.plot(t, z_wlc,  color='k')
        plt.title(f'z_bead = {z_bead_string} (um)')
        plt.ylim(bottom=0, top=np.max(z_wlc)*1.1)
        plt.hlines(-z_bead_opt, color='k', linestyle='--', xmin=t.min(), xmax=t.max())
        plt.tight_layout()
        plt.show()
    return z_bead_opt, z_bead_error

def get_lDNA(filename, label, z_bead):
    data = tio.hdf_data(filename, label)
    z = data.traces['Z (um)']
    f = data.traces['F (pN)']
    t = data.traces['Time (s)']


    z -= data.calculate_drift()
    z -= z_bead

    l = data.get_parameter('L (bp)')
    zWLC = fs.zWLC(f, l)

    mask = f > 5  # only consider forces above 5 pN

    lDNA = l* z/ zWLC

    levels = [5550, 4885,5686,5408,5348,5320,5047]


    h, bins = np.histogram(median_filter(lDNA[mask],10), bins=500, range=(0, l*1.1))
    bins = (bins[1:]+ bins[:-1])/2
    plt.plot(bins, h)
    plt.vlines(levels, 0, np.max(h), color='k', linestyle='--')
    plt.xlabel('L_DNA (bp)')
    plt.ylabel('Counts')
    plt.title(f'L = {l} (bp)')
    plt.show()



    plt.plot(f, lDNA, '.', color='blue', markersize=5, alpha=1)
    plt.plot(f, median_filter(lDNA,10), color='k', linewidth=2)
    plt.hlines(levels, 0, np.max(f), color='k', linestyle='--')
    plt.xlabel('F (pN)')
    plt.ylabel('L_DNA (bp)')
    plt.show()

    plt.plot(z, f, '.', color='blue', markersize=5, alpha=1)
    plt.plot(zWLC, f, color='k')
    for level in levels:
        plt.plot(fs.zWLC(f, level), f, color='k')
    plt.xlabel('Z (um)')
    plt.ylabel('F (pN)')
    plt.title(f'L = {l} (bp)')
    plt.show()


    return l, None

if __name__ == "__main__":
    label = '8'
    filename = Path("data/data_008.hdf")
    z_bead, _ = get_zbead(filename, label, show=False)

    get_lDNA(filename, label, z_bead)