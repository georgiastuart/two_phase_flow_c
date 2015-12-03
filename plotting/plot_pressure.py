import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    dim = 8

    pressure = np.reshape(np.fromfile('output/pressure.dat', dtype = np.float64), (dim, dim))

    title = "Pressure Field %d x %d" % (dim, dim)

    print pressure

    plt.pcolormesh(pressure)
    plt.title(title)
    plt.axis([0, dim, 0, dim])
    plt.gca().invert_yaxis()
    plt.show()
