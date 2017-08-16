#!/usr/bin/python3

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# over species
for i_field_idx in range(1,4):

    print('field = ',i_field_idx)

    # over timestep
    for t in range(0,500,100):

        print('time = ',t)

        basefname = 'u_' + str(i_field_idx) + '_' + str(t)
        datafname = '../reaction_system/' + basefname + '.csv'

        # x,y are spatial coords and z is concentration
        b = np.genfromtxt(datafname, skip_header=1, delimiter=',')
        x = b[:,1]
        y = b[:,2]
        z = b[:,0]
        
        # use this method to set aspect ratio because axis equal does not work with set_xbound, set_ybound (bug rep found on github)
        scale = 5
        dims = (scale*(max(x)-min(x)), scale*(max(y)-min(y)))
        f = plt.figure(figsize=dims)
        ax1 = f.add_subplot(111)
        
        ax1.tripcolor(x,y,z,cmap='jet')

        # add black circle for cylinder in the middle of flow
        circ = plt.Circle((0.20, 0.20), 0.05, color='k')
        ax1.add_artist(circ)
        
        ax1.set_xbound([0, max(x)])
        ax1.set_ybound([0, max(y)])
        plt.xlabel(r"$x$")
        plt.ylabel(r"$y$")
        
        imgfname = '../img/' + basefname + '.pdf'
        f.savefig(imgfname, bbox_inches='tight')
