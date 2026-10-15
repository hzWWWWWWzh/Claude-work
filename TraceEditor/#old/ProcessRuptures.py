import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from sklearn.linear_model import RANSACRegressor
from sklearn.linear_model import HuberRegressor
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
import os
from scipy.optimize import curve_fit


def calc_scatterplot(filename, treshold=0.5):
    xlabel = 'ln[(N$_0$-N+1/N) dF/dt (pN/s)]'
    ylabel = 'Force (pN)'
    if 'Telo' in filename:
        colors = [[1, 0, 0],
                  [0.286, 0.631, 0.243],
                  [0.243, 0.631, 0.392],
                  [0.243, 0.608, 0.631],
                  [1, 0, 0]]
    else:
        colors = [[0.631, 0.243, 0.592],
                  [0.243, 0.392, 0.631],
                  [0.369, 0.243, 0.631],
                  [0.431, 0.243, 0.631],
                  [0.631, 0.243, 0.592]]

    concentrations = [0, 1, 10, 100, 0]
    plot_x = [0, 1, 2, 3, 0.2]

    df = pd.read_excel(filename)
    dat = []
    for i, (x_col, f_rupt_col) in enumerate(np.reshape(np.asarray(df.columns), (-1, 2))):
        x = df[x_col].dropna().values
        y = df[f_rupt_col].dropna().values

        xmin, xmax = -4, 5
        ymin, ymax = 0, 45

        # Peform the kernel density estimate
        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([x, y])
        kernel = st.gaussian_kde(values)
        f = np.reshape(kernel(positions).T, xx.shape)

        fig = plt.figure(figsize=[4, 4])
        ax = fig.gca()
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        # Calculate levels
        n_levels = 5
        selected_levels = np.max(f)*np.linspace(treshold, 1, n_levels)

        # Contourf plot
        cmap = LinearSegmentedColormap.from_list(
            "white_to_rgb", [np.array([1, 1, 1]), np.asarray(colors[i])], N=256)
        cfset = ax.contourf(xx, yy, f, cmap=cmap, levels=selected_levels)
        cset = ax.contour(xx, yy, f, colors='k', levels=selected_levels)

        # Create a boolean array to check if points lie inside the top 3 contour levels
        # get paths of top 3 contour levels
        paths = [cset.collections[i].get_paths() for i in range(3)]

        selection = np.array([
            any(path.contains_point((xi, yi))
                for single_path_list in paths for path in single_path_list)
            for xi, yi in zip(x, y)])

        plt.scatter(x, y, s=0.1, color='k', alpha=0.5)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        n_total = len(x)
        x = x[selection]
        y = y[selection]

        p_lin, pcov = np.polyfit(x, y, 1, cov=True)
        e_lin = np.sqrt(np.diag(pcov))
        slope = p_lin[0]
        intercept = p_lin[1]

        line_X = np.linspace(xmin, xmax, 100)
        plt.plot(line_X,  intercept + slope*line_X,
                 color='k', linestyle='--', linewidth=2)
        plt.title(
            f'{concentrations[i]} nM TRF2$^\Delta$$ ^N$ (N={len(x)}/{n_total}) \n')

        folder_out = f'threshold{treshold}'
        if not os.path.exists(folder_out):
            os.makedirs(folder_out)

        plt.savefig(
            rf'{folder_out}/{filename.split(".")[0]}_{x_col}.jpg', bbox_inches='tight', dpi=600)
        plt.close()

        def f_rupture(x, d, k):
            kT = 4.114
            return (kT / d) * (x - np.log(kT * k / d))

        def g_rupture(x, d, G):
            kT = 4.114
            k0 = 1e7
            k = k0 * np.exp(-G)
            return (kT / d) * (x - np.log(kT * k / d))

        p_rupt, pcov = curve_fit(f_rupture, x, y, p0=[1, 1e-3])
        e_rupt = np.sqrt(np.diag(pcov))

        g_rupt, pcov = curve_fit(g_rupture, x, y, p0=[1, 25])

        # k0 = 1e9
        # G_barrier = -np.log(p_rupt[1] / k0)
        G_barrier = g_rupt[1]
        e_G_barrier = np.sqrt(np.diag(pcov))[1]

        dat.append([concentrations[i],
                   p_rupt[0], e_rupt[0],
                   p_rupt[1], e_rupt[1],
                   G_barrier, e_G_barrier,
                   p_lin[0], e_lin[0],
                   p_lin[1], e_lin[1],
                   plot_x[i], colors[i], x_col])

    df = pd.DataFrame(
        dat, columns=['c (nM)', 'd (nm)', 'ste[d (nm)]', 'kd (s-1)', 'ste[kd (s-1)]', 'G_barrier (kT)', 'ste[G_barrier (kT)]', 'slope', 'ste[slope]', 'intercept', 'ste[intercept]', 'x', 'colors', 'column'])
    df['treshold'] = treshold
    df['filename'] = filename

    return df, folder_out


def plot_data(df, name, filename=None, lim=None):

    fig, ax = plt.subplots(figsize=(4, 4))

    for i, row in df.iterrows():
        fillstyle = 'none' if 'trol' in row['column'] else 'full'
        try:
            ax.errorbar(row['x'], row[name], yerr=row[f'ste[{name}]'], color=row['colors'],
                        fmt='o', capsize=8, fillstyle=fillstyle)
        except KeyError:
            ax.plot(row['x'], row[name], 'o', color=row['colors'],
                    fillstyle=fillstyle)
    add_ticks = ['control' not in colx for colx in df['column']]
    ax.set_xticks(df[add_ticks]['x'])
    ax.set_xticklabels(df[add_ticks]['c (nM)'])
    if name == 'kd (s-1)':
        name = 'k$_0$ (s$^{-1}$)'
    elif name == 'G_barrier (kT)':
        name = 'G$^{\dag}$ ($k_B$T)'

    ax.set_ylabel(name)
    ax.set_xlabel(f'[TRF2$^\Delta$$ ^N$] (nM)')
    if lim is not None:
        ax.set_ylim(lim)

    if filename is not None:
        fig.tight_layout()
        fig.savefig(filename, bbox_inches='tight', dpi=600)
    else:
        plt.show()
    return


filenames = [r'Ruptures-Telo.xlsx', r'Ruptures-601.xlsx']

# for treshold in [0, 0.25, 0.5, 0.75]:
for treshold in [0.5]:
    df = pd.DataFrame()
    for filename in filenames[::-1]:
        df_new, folder_out = calc_scatterplot(filename, treshold)
        df = pd.concat([df, df_new], ignore_index=True)

    print(df)
    df.to_excel(rf'{folder_out}/Summary.xlsx', index=False)
    print(rf'Saved to {folder_out}/Summary.xlsx')

    for name, lim in {'d (nm)': [0, 2.5], 'kd (s-1)': [0, 0.02], 'G_barrier (kT)': [24, 28]}.items():
        plot_data(df, name, filename=rf'{folder_out}/{name}.jpg', lim=lim)
