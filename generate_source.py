import numpy as np

dim = 128
source = np.zeros([dim, dim])
source[dim - 1][0] = 10**(-7)
source[0][dim - 1] = -10**(-7)

np.savetxt('input/src_field_large.txt', source, delimiter='\n')

perm = np.zeros([dim, dim])
np.savetxt('input/perm_field_large.txt', perm, delimiter='\n')
