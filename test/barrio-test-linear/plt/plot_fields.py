#!/usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

subdir = 'pnl_d/'

for t in range(0,2001,50):
    fname = subdir + 't' + str(t) + '_f0'
    datafname = 'dat/' + fname + '.dat'
    field = np.loadtxt(datafname)
    plt.imshow(field, cmap='jet', interpolation='nearest')
    plt.colorbar()
    imgfname = 'img/' + fname + '.pdf'
    plt.savefig(imgfname, bbox_inches='tight')
    plt.clf()

# plt.show()
