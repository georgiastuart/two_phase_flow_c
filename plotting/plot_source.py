import matplotlib.pyplot as plt
import numpy as np

source = np.loadtxt('output/source_check.txt', delimiter='\n')
sat = np.arange(.21, .841, 0.001)

print source.shape
print sat.shape
plt.plot(sat, source)
plt.show()
