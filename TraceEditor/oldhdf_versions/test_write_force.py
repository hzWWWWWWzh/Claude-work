
import numpy as np
from matplotlib import pyplot as plt

from icecream import ic

import TraceIO as tio
from ForceSpectroscopy import MagnetForce

def test_write_force(filename):
    data = tio.hdf_data(filename)
    ic(data.list_parameters())

    data.set_parameter('Fmax (pN)', 100, type='global')
    data.set_parameter('Fmin (pN)', 0.02, type='global')
    data.set_parameter('L1 (mm)', 1.4, type='global')
    data.set_parameter('L2 (mm)', 0.8, type='global')
    data.set_parameter('alpha (-)', 0.7, type='global')
    

    force = MagnetForce(data.traces['Shift (mm)'].values,
                        fmax=data.get_parameter('Fmax (pN)'),
                        f0=data.get_parameter('Fmin (pN)'),
                        L1=data.get_parameter('L1 (mm)'),
                        L2=data.get_parameter('L2 (mm)'),
                        alpha=data.get_parameter('alpha (-)'))

    data.traces['F (pN)'] = force
    data.save(traces=True, parameters=True)
    plt.show()

def plot_shift_force(filename):
    data = tio.hdf_data(filename)
    shift = data.traces['Shift (mm)'].values
    force = data.traces['F (pN)'].values

    for p in ['Fmax (pN)', 'Fmin (pN)', 'L1 (mm)', 'L2 (mm)', 'alpha (-)']:
        ic(p, data.get_parameter(p))


    plt.plot(shift, force)
    plt.xlabel('Shift (mm)')
    plt.ylabel('F (pN)')
    plt.yscale('log')
    plt.xticks(np.arange(0, 15, 1))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    filename = r"data\data_021.hdf"
    test_write_force(filename)
    plot_shift_force(filename)