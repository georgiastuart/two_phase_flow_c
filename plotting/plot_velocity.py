import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    dim = 8

    vel_y = np.reshape(np.fromfile('output/velocity_y.dat', dtype = np.float64), (dim, dim))
    vel_x = np.reshape(np.fromfile('output/velocity_x.dat', dtype = np.float64), (dim, dim))

    title = "Velocity Field %d x %d" % (dim, dim)

    plt.quiver(vel_x, vel_y)
    plt.title(title)
    plt.axis([0, dim - 1, 0, dim - 1])
    plt.gca().invert_yaxis()
    plt.show()
