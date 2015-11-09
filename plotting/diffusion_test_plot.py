from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

if __name__ == '__main__':
    config = ConfigParser.ConfigParser()
    config.read('input/config.ini')
    dim = float(config.get('dimensions', 'xdim'))
    h = 1/dim

    print h

    temp = np.reshape(np.fromfile('output/saturation.dat', dtype = np.float64), (dim, dim))
    temp_1D = temp[1][:]
    x = np.arange(0, 1, h)
    analytical = np.exp(-np.pi**2 * h) * np.sin(np.pi * x) + 0.5 * np.exp(-9 * np.pi**2 * h) * np.sin(3 * np.pi * x)
    ic = np.sin(np.pi * x) + 0.5 * np.sin(3 * np.pi * x)

    plt.plot(x, temp_1D, label = 'Numerical Solution')
    plt.plot(x, analytical, label = 'Analytical Solution')
    plt.plot(x, ic, label = 'Initial Condition')
    plt.legend()

    plt.show()
