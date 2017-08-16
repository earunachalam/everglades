#!/usr/bin/python

# 1d.py
# attempt to reproduce Fig.1 of Maini, et al. Interface Focus (2012)
# start from uniform state and attempt to reach the spatially periodic steady state


from __future__ import print_function
from fenics import *
from mshr import *
import numpy as np

import sys

# final time
T = 0.5

# time step size
dt = 0.01

# number of timesteps
nsteps = int(T / dt)
print nsteps, " steps total. Starting..."

# parameters for Barrio model
alpha = 0.899
beta = -0.91
gamma = -0.899
delta = 2.0
r1 = 0.02
r2 = 0.2
d = 0.536

class InitialCondition(Expression):
    def eval(self, value, x):
        value[0] = 1#np.random.rand()
        value[1] = 1#np.random.rand()
    def value_shape(self):
        return (2,)

# Create VTK files for visualization output
basedir = 'dat/fish-pigmentation/'
vtkfile_u_1 = File(basedir + 'u_1.pvd')
vtkfile_u_2 = File(basedir + 'u_2.pvd')

# Define function space for system of concentrations
mesh = UnitSquareMesh(100, 100)
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2 = TestFunctions(V)

# Define functions for concentrations
u = Function(V)
u_n = interpolate(InitialCondition(degree=2), V)

# Split system functions to access components
u_1, u_2 = split(u)
u_n1, u_n2 = split(u_n)

# plot(u_n2)
# interactive()
# sys.exit()

# Define source terms
f_1 = Constant(0)
f_2 = Constant(0)

# Define expressions used in variational forms
k = Constant(dt)

alpha = Constant(alpha)
beta = Constant(beta)
gamma = Constant(gamma)
delta = Constant(delta)
r1 = Constant(r1)
r2 = Constant(r2)
d = Constant(d)

# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx + delta*d*dot(grad(u_1), grad(v_1))*dx + (alpha*u_1*(1.0-r1*u_2*u_2) + u_2*(1.0-r2*u_1))*v_1*dx \
  + ((u_2 - u_n2) / k)*v_2*dx + delta*dot(grad(u_2), grad(v_2))*dx + (beta*u_2*(1.0+r1*u_1*u_2/beta) + u_1*(gamma+r2*u_2))*v_2*dx \
  - f_1*v_1*dx - f_2*v_2*dx

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

# Time-stepping
t = 0
for n in range(nsteps):

    # Update current time
    t += dt

    # Solve variational problem for time step
    solve(F == 0, u)

    # Save solution to file (VTK)
    _u_1, _u_2 = u.split()
    vtkfile_u_1 << (_u_1, t)
    vtkfile_u_2 << (_u_2, t)

    # Update previous solution
    u_n.assign(u)

    # Update progress bar
    progress.update(t / T)

