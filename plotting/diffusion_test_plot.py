from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

if __name__ == '__main__':
    config = ConfigParser.ConfigParser()
    config.read('input/config.ini')
    dim = float(config.get('dimensions', 'xdim'))
    h = 1/dim
    ts = h / 2

    print h

    temp = np.reshape(np.fromfile('output/saturation_coarse.dat', dtype = np.float64), (dim, dim))
    temp_1D_coarse = temp[1][:]
    temp = np.reshape(np.fromfile('output/saturation.dat', dtype = np.float64), (dim * 2, dim * 2))
    temp_1D_fine = temp[1][:]
    x = np.arange(0, 1, h)
    x_fine = np.arange(0, 1, h/2)
    analytical = np.exp(-np.pi**2 * ts) * np.sin(np.pi * x) + 0.5 * np.exp(-9 * np.pi**2 * ts) * np.sin(3 * np.pi * x)
    ic = np.sin(np.pi * x) + 0.5 * np.sin(3 * np.pi * x)

    plt.plot(x, temp_1D_coarse, label = 'Numerical Solution (Coarse)')
    plt.plot(x_fine, temp_1D_fine, label = 'Numerical Solution (Fine)')
    plt.plot(x, analytical, label = 'Analytical Solution')
    plt.plot(x, ic, label = 'Initial Condition')
    plt.legend()

    plt.show()
