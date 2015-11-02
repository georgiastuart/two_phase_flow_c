import numpy as np

if __name__ == '__main__':
    pressure = np.fromfile('output/pressure.dat', dtype = np.float64)

    pressure_reshape = pressure.reshape((64, 64))
    print pressure_reshape

    np.savetxt('pressure_plaintxt.txt', pressure_reshape, delimiter = '\t', newline = '\n')
