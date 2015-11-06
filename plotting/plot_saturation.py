import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    dim = 8

    pressure = np.reshape(np.fromfile('output/saturation.dat', dtype = np.float64), (dim, dim))

    title = "Saturation Field %d x %d" % (dim, dim)

    plt.imshow(pressure)
    plt.title(title)
    plt.colorbar()
    plt.axis([0, dim - 1, 0, dim - 1])
    plt.gca().invert_yaxis()
    plt.show()
