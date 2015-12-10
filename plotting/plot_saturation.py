import numpy as np
import matplotlib.pyplot as plt
import ConfigParser as cfg

if __name__ == '__main__':
    parser = cfg.ConfigParser()
    parser.read('input/config.ini')
    dim = parser.getint('dimensions', 'xdim')

    print dim

    saturation = np.reshape(np.fromfile('output/saturation.dat', dtype = np.float64), (dim, dim))

    title = "Saturation Field %d x %d" % (dim, dim)

    plt.pcolormesh(saturation)
    plt.title(title)
    plt.colorbar()
    plt.axis([0, dim, 0, dim])
    plt.gca().invert_yaxis()
    plt.show()

    plt.contour(saturation)
    plt.title(title)
    plt.colorbar()
    plt.axis([0, dim, 0, dim])
    plt.gca().invert_yaxis()
    plt.show()
