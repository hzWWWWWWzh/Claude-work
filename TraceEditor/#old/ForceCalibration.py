import numpy as np
import matplotlib.pyplot as plt
from scipy import signal, optimize
import TraceIO as tio
import pandas as pd

kT = 4.1e-21


def Brownian_dynamics(tmax, eta=1e-3, R=1.4e-6, k=1e-7, fs=60, W=0.016, sigma=10, show=False):
    """
    Create a timetrace according to Wong2006OptExpress equation 8
    :param tmax: duration (s)
    :param eta: Viscosity (N/m2s)
    :param R: Bead radius (m)
    :param k: Stiffness of bead trap (N/m)
    :param fs: Sampling frequency (Hz)
    :param W: Integration time (s)
    :param sigma: tracking accuracy (nm)
    :param show: Plot the curve?
    :return: time (s), position (um)
    """
    n_samples = int(tmax * fs)
    n_subsamples = 1000
    dt = 1 / fs
    gamma = 6 * np.pi * eta * R
    D = kT / gamma
    tau = gamma / k
    varx = 2 * D / (fs * n_subsamples)

    x_sub = np.zeros(n_subsamples)
    x = 0
    x_t = []
    for i in range(n_samples):
        for i, dx in enumerate(np.random.normal(0, np.sqrt(varx), n_subsamples)):
            x -= x * (dt / n_subsamples) / tau + dx
            x_sub[i] = x
        x_t.append(np.mean(x_sub[:-int(W * fs * n_subsamples)]))
    x_t = np.asarray(x_t) + np.random.normal(0, sigma * 1e-9, n_samples)
    tmax = np.linspace(0, tmax, n_samples)
    if show:
        plt.plot(tmax, x_t * 1e6, color='red')
        tio.format_plot('t (s)', 'x (um)', yrange=[-1, 1])

    return tmax, x_t * 1e6


def fit_powerspectrum(t, x, show=False, eta=1e-3, R=5e-7, W=0.016):
    """
    Fit powerspectrum according to Wong2006OptExpress. Weigh spectrum proportional to frequency
    to reduce impact of poorly measured low frequencies.
    :param t: time (s)
    :param x: position (um)
    :param show: Show the plot?
    :param eta: Viscosity (N/m2s)
    :param R: Bead radius (m)
    :param W: Integration time (s)
    :return: stiffness (N/m)
    """
    fs = 1 / np.median(np.diff(t))
    R = 1.4e-6

    def powerspectrum(f, k=1e-7, sigma=10):
        """
        Fit function for powerspectrum, including aliasing and blurring.
        Add more parameters in function io if you want to fit them,
        otherwise parameters are taken from global space.
        :param f: array of frequencies (Hz)
        :param k: stiffness (N/m)
        :param sigma: tracking accuracy (nm)
        :return: power density (um2/Hz)
        """
        gamma = 6 * np.pi * eta * R
        spectrum = np.zeros_like(f)
        for n in [-1, 0, 1]:
            spectrum += 4 * gamma * kT / ((gamma ** 2) * (((n * fs + f) * 2 * np.pi) ** 2) + k ** 2) \
                        * np.sinc((n * fs + f) * W) ** 2
        spectrum += (sigma * 1e-9) ** 2
        return spectrum * 1e12

    freqs, Pxx = signal.welch(x, fs, 'hanning', scaling='density')
    names = ['k (N/m)', 'sigma (nm)']
    p0 = [1e-6, 10]
    params, cov = optimize.curve_fit(powerspectrum, freqs, Pxx, sigma=1 / freqs, p0=p0)

    for i, p in enumerate(params):
        print(tio.round_significance(np.abs(p), np.sqrt(cov[i, i]), par_name=names[i]))

    fit = powerspectrum(freqs, *params)
    if show:
        plt.scatter(freqs, Pxx, color='red')
        plt.plot(freqs, fit, color='black')
        tio.format_plot('f (Hz)', 'x2 (um2/Hz)', ylog=True, xlog=True, xrange=[1e-1, fs / 2], yrange=[1e-5, 1e0])
    return params


if __name__ == "__main__":
    filename = r'c:\tmp\test.dat'
    do_brownian_dynamics = True
    if do_brownian_dynamics:
        t, x = Brownian_dynamics(240, fs=60, k=3.2e-6, sigma=10, show=True)
        df = pd.DataFrame({'t (s)': t, 'x (um)': x})
        df.to_csv(filename)
    else:
        df = pd.read_csv(filename, index_col=[0])
        t = df['t (s)'].to_numpy()
        x = df['x (um)'].to_numpy()

    fit_powerspectrum(t, x, show=True)
