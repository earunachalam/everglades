#!/usr/bin/python3

from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack as fp
import scipy.stats as sc
import sys

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

ref_dist = 0
for i in ['a','b','c','d']:
    psd2D = np.zeros([100,100])
    for ts in range(1900,2000,1):
        data_x = np.loadtxt('dat/pnl_' + i + '/t' + str(ts) + '_f0.dat')
        data_k = fp.fftshift(fp.fft2(data_x))
        psd2D += np.abs(data_k)**2
        freq = np.fft.fftfreq(data_x.shape[-1])
    
    if i == 'a':
        # make reference state truly homogeneous
        ref_dist = np.copy(psd2D)
        # ref_dist.fill(1.)           # does not matter what value is used to fill
        ref_dist.fill(np.finfo(float).eps)
        print(psd2D)
    
    print(i,"/ homog KL div = ", sc.entropy(psd2D.flatten(),ref_dist.flatten()))
    
    n = np.min(freq)
    x = np.max(freq)
    ct = [n,x,n,x]
    
    psd2D[psd2D == np.max(psd2D)] = np.min(psd2D)
    plt.imshow(psd2D, cmap='gnuplot', interpolation='nearest', extent=ct)
    plt.colorbar()
    ax = plt.gca()
    ax.set_xlim(-0.15,0.15)
    ax.set_ylim(-0.15,0.15)
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')

    plt.savefig("fft_" + i + '.pdf', bbox_inches='tight')
    plt.close()
