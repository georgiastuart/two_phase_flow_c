import numpy as np

source = np.zeros([64, 64])
source[63][0] = 10**(-7)
source[0][63] = -10**(-7)

np.savetxt('src_field.txt', source, delimiter='\n')
