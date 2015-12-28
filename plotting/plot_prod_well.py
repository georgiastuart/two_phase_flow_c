import numpy as np
import matplotlib.pyplot as plt
import ConfigParser as cfg

if __name__ == '__main__':
    parser = cfg.ConfigParser()
    parser.read('input/config.ini')
    ts = parser.getint('dimensions', 'time_steps')

    data = np.fromfile('output/prod_well.dat', dtype=np.float64);

    print data.shape
    print data;
