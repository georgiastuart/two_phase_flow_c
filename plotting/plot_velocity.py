import numpy as np
import matplotlib.pyplot as plt
import ConfigParser as cfg

if __name__ == '__main__':
    parser = cfg.ConfigParser()
    parser.read('input/config.ini')
    dim = parser.getint('dimensions', 'xdim')

    vel_y = np.reshape(np.fromfile('output/velocity_y.dat', dtype = np.float64), (dim, dim))
    vel_x = np.reshape(np.fromfile('output/velocity_x.dat', dtype = np.float64), (dim, dim))

    title = "Velocity Field %d x %d" % (dim, dim)

    plt.quiver(vel_x, vel_y)
    plt.title(title)
    plt.axis([-1, dim, -1, dim])
    plt.gca().invert_yaxis()
    plt.show()
