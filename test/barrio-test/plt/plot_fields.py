#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np

field = np.loadtxt('/home/easun/Research/everglades/test/barrio-test/dat/pnl_c/t999_f0.dat')
plt.imshow(field, cmap='jet', interpolation='nearest')

plt.show()
