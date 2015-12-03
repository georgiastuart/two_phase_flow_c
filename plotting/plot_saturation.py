import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    dim = 8

    saturation = np.reshape(np.fromfile('output/saturation.dat', dtype = np.float64), (dim, dim))

    title = "Saturation Field %d x %d" % (dim, dim)

    plt.pcolormesh(saturation)
    plt.title(title)
    plt.colorbar()
    plt.axis([0, dim, 0, dim])
    plt.gca().invert_yaxis()
    plt.show()
