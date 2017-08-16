#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys

try:
    fn_settings = sys.argv[1]
    fn_output = sys.argv[2]
except:
    print("Syntax: initial_condition.py filename_settings filename_output")

d = {}
with open(fn_settings) as f:
    for line in f:
       (key, val) = line.split()
       d[key] = val

Nx1         = int(d['Nx1'])
Nx2         = int(d['Nx2'])
shape       = (Nx1,Nx2)
initcond    =  float(d['init_base_amp'])*np.ones(shape, dtype=np.float)
initcond    += float(d['init_rand_amp'])*np.random.rand(Nx1, Nx2)

# print(initcond)
# plt.imshow(initcond, cmap='jet', interpolation='nearest')
# plt.show()
