import ForceSpectroscopy as fs
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import product
from tqdm import tqdm

KT = 4.1  # pN nm
OMEGA_0 = 2 * np.pi / 3.6  # rad/nm


def free_energy(force, sigma, state="twisted"):
    if state == "twisted":
        c_twisted = 100  # nm
        a_twisted = 50  # nm
        c_twisted *= 1 - (c_twisted / (4 * a_twisted)) * np.sqrt(
            1 / (a_twisted * force)
        )
        g_state = -(force / KT - np.sqrt((force / KT) / a_twisted))
        g_state += 0.5 * c_twisted * sigma**2 * OMEGA_0**2
    elif state == "plectonemic":
        c_plectonemic = 24  # nm
        g_state = 0.5 * c_plectonemic * sigma**2 * OMEGA_0**2
    elif state == "melted":
        c_melted = 1.75  # nm
        a_melted = 4  # nm
        sigma_0 = -1
        epsilon = 2  # kT/nm
        g_state = -((force / KT) - np.sqrt((force / KT) / a_melted)) + epsilon
        g_state += 0.5 * c_melted * (sigma - sigma_0) ** 2 * OMEGA_0**2
    return g_state


def plot_free_energy(force):
    sigma = np.linspace(-1.5, 1, 1000)
    for state in ["twisted", "plectonemic", "melted"]:
        g = free_energy(force, sigma, state)
        plt.plot(sigma, g, label=state)
    plt.legend()
    plt.ylim(-0.9, 7)
    plt.xlim(-1.5, 0.5)
    plt.xlabel("sigma")
    plt.ylabel("Free energy (kT/nm)")
    plt.show()


def plot_force_extension(L=1000):
    f = np.logspace(-1.4, 2, 1000)
    plt.plot(1000 * fs.zWLC(f, L), f, label="dsDNA")
    plt.plot(fs.zFJC(f / 2, b=1.55, L=1.6 * L / 3), f, label="ssDNA")
    plt.xlabel("z (nm)")
    plt.ylabel("Force (pN)")
    plt.legend()
    plt.show()


def calc_extension(L=3000, force=1.0, sigma=-0.05, steps=64):
    # make 4d array of all possible combinations
    n_t, n_p, s_t, s_p = np.meshgrid(
        np.linspace(0, L, steps),
        np.linspace(0, L, steps),
        np.linspace(-0.1, 0.1, steps),
        np.linspace(-0.1, 0.1, steps),
        indexing="ij",
    )

    # calculate the number of melted base pairs, and discard invalid combinations
    n_m = L - n_t - n_p
    valid = n_m > 0
    n_t, n_p, n_m, s_t, s_p = n_t[valid], n_p[valid], n_m[valid], s_t[valid], s_p[valid]

    if isinstance(sigma, (int, float)) and isinstance(force, np.ndarray):
        sigma = np.full_like(force, sigma)
    elif isinstance(sigma, np.ndarray) and isinstance(force, (int, float)):
        force = np.full_like(sigma, force)

    extension = np.zeros_like(force)
    for i, (f, s) in enumerate(tqdm(zip(force, sigma), total=len(force))):
        # calculate the twist of melted base pairs
        s_m = (L * s - s_t * n_t - s_p * n_p) / n_m
        valid = (s_m >= -0.1) & (s_m <= 0)

        z = fs.zWLC(f, n_t[valid]) * 1000
        +fs.zFJC(f / 2, b=1.55, L=1.6 * n_m[valid] / 3)

        energy = free_energy(f, s_t[valid], "twisted") * n_t[valid]
        +free_energy(f, s_p[valid], "plectonemic") * n_p[valid]
        +free_energy(f, s_m[valid], "melted") * n_m[valid]
        -f * z

        exp_energy = np.exp(-(energy - np.min(energy)))
        extension[i] = np.sum(z * exp_energy) / np.sum(exp_energy)

    return extension


if __name__ == "__main__":
    force = np.logspace(-1.4, 1.8, 100)
    # plot_free_energy(1)
    # plot_force_extension(L=3000)

    sigma = np.linspace(-0.1, 0.1, 50)
    z = calc_extension(sigma=sigma, force=2.95, L=1000, steps=100)

    plt.plot(sigma, z, "o-")
    plt.show()
