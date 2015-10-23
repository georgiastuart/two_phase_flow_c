import numpy as np

dim = 8
source = np.zeros([dim, dim])
# source[dim - 1][0] = 10**(-7)
# source[0][dim - 1] = -10**(-7)

np.savetxt('perm_field_small.txt', source, delimiter='\n')
