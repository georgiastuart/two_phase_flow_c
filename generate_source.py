import numpy as np

dim = 64
suffix = "64"
source = np.zeros([dim, dim])
source[dim - 1][0] = 10**(-7)
source[0][dim - 1] = -10**(-7)

src_name = 'input/src_field_%s.txt' % (suffix)
np.savetxt(src_name, source, delimiter='\n')

perm_name = 'input/perm_field_%s.txt' % (suffix)
perm = np.zeros([dim, dim])
np.savetxt(perm_name, perm, delimiter='\n')

sat_name = 'input/sat_field_%s.txt' % (suffix)
sat = 0.21 * np.ones([dim, dim])
sat[dim - 1][0] = 0.84
np.savetxt(sat_name, sat, delimiter='\n')
