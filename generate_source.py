import numpy as np
import argparse as ap

parser = ap.ArgumentParser(description = 'Generate testing fields')
parser.add_argument('dimension', metavar='D', type=int, help="Dimension of the grid")
parser.add_argument('src_dim', metavar='S', type=int, help='Dimension of the source')

args = parser.parse_args()
dim = args.dimension
source = np.zeros([dim, dim])
source_size = args.src_dim

for i in xrange(source_size):
    for j in xrange(source_size):
        source[dim - 1 - i][j] = 10**(-7)
        source[i][dim - 1 - j] = -10**(-7)

src_name = 'input/src_field_%d.txt' % (dim)
np.savetxt(src_name, source, delimiter='\n')

perm_name = 'input/perm_field_%d.txt' % (dim)
perm = np.zeros([dim, dim])
np.savetxt(perm_name, perm, delimiter='\n')

sat_name = 'input/sat_field_%d.txt' % (dim)
sat = 0.21 * np.ones([dim, dim])
# for i in xrange(source_size):
#     for j in xrange(source_size):
#         sat[dim - 1 - i][j] = 0.84
np.savetxt(sat_name, sat, delimiter='\n')
