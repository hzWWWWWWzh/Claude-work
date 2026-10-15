import ForceSpectroscopy as fs
import TraceIO as tio
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy import sparse
from time import time

L_bp = 30


def zDNA_noneq2(time, f_pN, show=True):
    b_nm = 1.55
    Ebp_kT = 1.57
    n_strands = 2

    # calculate extension and free energy as a function of force for each state
    extension = [fs.zWLC(f_pN),
                 1.15 * fs.zFJC(f_pN / n_strands, b=b_nm, EFJC=True)]
    # JvN: Why do I need 1.15 fudgefactor ?????
    free_energy = [fs.gWLC(f_pN),
                   n_strands * fs.gEFJC(f_pN / n_strands, b=1.5) + Ebp_kT]

    gmatrixes, zstates = calc_gmatrixes(f_pN, extension, free_energy)

    px0 = np.exp(-np.diagonal(gmatrixes[0]))
    px0 /= np.sum(px0)
    px = [px0]
    for dt, gmatrix in zip(np.diff(time), gmatrixes[:-1]):
        pmatrix = calc_pmatrix(calc_tmatrix(gmatrix), dt)
        px.append(pmatrix.T.dot(px[-1]))

    if show:
        plt.plot(f_pN, px)
        tio.format_plot('F (pN)', 'P(State)', yrange=[-0.1, 1.1])

    #   plt.plot(f_pN, zstates)
    #   plt.plot(f_pN, np.sum(px*zstates, axis=1))
    #   tio.format_plot('F (pN)', 'z (nm)', yrange = [-0.1, 1.1])

        plt.show()
    return


def zDNA_noneq(f_pN, show=True):
    b_nm = 1.55
    Ebp_kT = 1.57
    n_strands = 2

    zds = fs.zWLC(f_pN)
    # JvN: Why do I need 1.15 fudgefactor ?????
    zss = 1.15*fs.zFJC(f_pN / n_strands, b=b_nm, EFJC=True)
    gds = fs.gWLC(f_pN)
    gss = n_strands * fs.gEFJC(f_pN / n_strands, b=1.5) + Ebp_kT
    zb, gb0 = fs.intersection(zds, gds, zss, gss)

    if show:  # plot potential energy landscape
        f = 30
        gds0 = np.min(gds - zds*f/fs.kT)
        zds0 = zds[np.argmin(gds - zds*f/fs.kT)]
        gss0 = np.min(gss - zss*f/fs.kT)
        zss0 = zss[np.argmin(gss - zss*f/fs.kT)]

        plt.plot(zds, gds - zds*f/fs.kT, color='blue')
        plt.plot(zss, gss - zss*f/fs.kT, color='red')
        plt.plot(zds0, gds0, 'o', color='blue', fillstyle='none', label='ds')
        plt.plot(zss0, gss0, 'o', color='red', fillstyle='none', label='ss')
        plt.plot(zb, gb0 - zb*f/fs.kT, 'o', color='black',
                 fillstyle='none', label='barier')
        plt.legend(loc='upper right')
        tio.format_plot('z (nm)', 'g (kT)', xrange=[
                        0, 0.75], aspect=1, scale_page=0.5, title=f'F = {f:.2f} (pN))')
        plt.show()

    gds0 = []
    gss0 = []
    gb = gb0 - zb*f_pN/fs.kT
    for f in f_pN:
        gds0.append(np.min(gds - zds*f/fs.kT))
        gss0.append(np.min(gss - zss*f/fs.kT))

    if show:  # Plot energies as a function of force
        plt.plot(f_pN, gds0, color='blue', label='ds')
        plt.plot(f_pN, gss0, color='red', label='ss')
        plt.plot(f_pN, gb, color='black', label='barier')
        plt.legend()
        tio.format_plot('F (pN)', 'g (kT)', aspect=1, scale_page=0.5)
        plt.show()


def calc_gmatrixes(force_pN, extension_nm, free_energy_kT):
    n_states = len(extension_nm)
    g_barier = np.zeros([n_states, n_states])
    z_barier = np.zeros_like(g_barier)
    for i, (z, g) in enumerate(zip(extension_nm, free_energy_kT)):
        for j, (z_j, g_j) in enumerate(zip(extension_nm, free_energy_kT)):
            if j < i:
                z_barier[i, j], g_barier[i, j] = fs.intersection(
                    [z, g], [z_j, g_j])
                z_barier[j, i] = z_barier[i, j]
                g_barier[j, i] = g_barier[i, j]
    gmatrices = []
    zstates = []
    for f in force_pN:
        g_matrix = g_barier - f * z_barier/fs.kT
        z_state = np.zeros(n_states)
        for i, (z, g) in enumerate(zip(extension_nm, free_energy_kT)):
            g_f = g - f * z / fs.kT
            p = np.exp(-g_f)
            p /= np.sum(p)
            g_matrix[i, i] = np.sum(p * g_f)
            z_state[i] = np.sum(p * z)
        gmatrices.append(g_matrix)
        zstates.append(z_state)
    return np.asarray(gmatrices), np.asarray(zstates)


def calc_timestep(gmatrix, p0=None, t_s=0.1, k_s=1e-5, single_molecule=False, show=False):
    """
    Compute a new probability distribution after dt seconds based on matrix of 
    free enegeries, bariers and privious probability distribution.
    Use sparse matrix manipulations to speed up large matrix computations.
    Steps:
    1) delta_g-matrix : delta_g = g_barier - g0
    2) t-matrix : t = ks * exp(- delta_g)
    3) Impose detailed balance: sum of transfer rate rows = 0
    4) p-matrix: p = matrixexponential(t * dt)
    5) new distribution of states: p(t) = T(p-matrix).p0

    For easy introduction see Jones2017MedDescMaking, DOI: 10.1177/0272989X17696997
    """
    if False:  # Convert to sparse g-matrix
        gmatrix = sparse.csc_matrix(gmatrix)
    if sparse.issparse(gmatrix):
        matrix = gmatrix
        matrix.data -= np.take(matrix.diagonal(), matrix.indices)
        matrix.data = k_s * np.exp(matrix.data)
        matrix -= sparse.dia_matrix((matrix.sum(axis=1).A1,
                                     [0]), np.shape(gmatrix))
        matrix = sparse.linalg.expm(matrix*t_s)
    else:
        matrix = np.zeros_like(gmatrix)
        for i, g in enumerate(gmatrix):
            matrix[i] = k_s * np.exp(g - g[i])
            matrix[i, i] = 0
            matrix[i, i] = -np.sum(matrix[i])
        matrix = linalg.expm(matrix * t_s)

    if p0 is None:
        p0 = np.zeros(len(gmatrix))
        p0[0] = 1

    if single_molecule:
        pt = np.zeros(len(p0))
        states = np.arange(len(p0), dtype=int)
        x = np.random.choice(states, p=matrix[np.argmax(p0)])
        pt[x] = 1
        return pt
    else:
        return matrix.T.dot(p0)


def print_matrix(matrix, title=None, precision=2):
    if title is not None:
        print(title)
    print(np.array_str(matrix, precision=precision, suppress_small=True))
    print()


if __name__ == "__main__":

    gmatrix = np.array([[0, 7, 9, 8, 9, 8, 9],
                        [7, 0, 9, 10, 7, 9, 9],
                        [9, 9, 2, 11, 11, 11, 11],
                        [8, 10, 11, 0, 8, 9, 10],
                        [9, 7, 11, 8, 1, 9, 10],
                        [8, 9, 11, 9, 9, -1, 8],
                        [9, 9, 11, 10, 10, 8, 0]], np.float32)
    gmatrix = np.array([[0, 10], [10, 1.5]])
    # gmatrix = np.array([[0, 10, 20], [10, 1.5, 10], [20, 10, 3]])

    # gmatrix= np.reshape(np.arange(9), (3,3))
    # gmatrix = np.zeros((10,10))
    # np.fill_diagonal(gmatrix,3)

    if True:  # test use of sparse matrix
        print_matrix(gmatrix, 'G-matrix:\n(diagonal = energy minima, non-diagonal are transition bariers)')
        time = np.arange(0, 10, 0.001)
        p0 = np.zeros(len(gmatrix))
        p0[-1] = 1
        p = [p0]
        show = True
        for dt in np.diff(time):
            p.append(calc_timestep(gmatrix, p[-1], t_s=dt, single_molecule=False))
            show = False
        plt.plot(time, p)
        tio.format_plot('time (s)', 'P(x)')
        plt.show()

    if False:  # Calculate non-equilibrium f-z curve
        shift = np.linspace(12, 0, 500)
        shift = np.append(shift, np.flip(shift))
        force = fs.MagnetForce(shift)
        time = np.linspace(0, 10, len(force))
        zDNA_noneq2(time, force, True)
