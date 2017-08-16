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
print(nsteps, " steps total. Starting...")
sys.stdout.flush()

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
        value[0] = 0.10*2.00*(np.random.rand()-0.50)
        value[1] = 0.10*2.00*(np.random.rand()-0.50)
    def value_shape(self):
        return (2,)

# Create VTK files for visualization output
basedir = 'dat/fish-pigmentation/'
vtkfile_u_1 = File(basedir + 'u_1.pvd')
vtkfile_u_2 = File(basedir + 'u_2.pvd')

Lx = 100
Ly = 100
Nx = 100
Ny = 100
lowerLeft = Point(0.0,0.0)
upperRight = Point(Lx,Ly)
mesh = RectangleMesh(lowerLeft, upperRight, Nx, Ny)

# class PeriodicDomain(SubDomain):

    # def inside(self, x, on_boundary):
        # return bool((near(x[0], 0) or near(x[1], 0)) and 
            # (not ((near(x[0], Lx) and near(x[1], 0)) or 
                  # (near(x[0], 0) and near(x[1], Ly)))) and on_boundary)

    # def map(self, x, y):
        # if near(x[0], Lx) and near(x[1], Ly):
            # y[0] = x[0] - Lx
            # y[1] = x[1] - Ly
        # elif near(x[0], Lx):
            # y[0] = x[0] - Lx
            # y[1] = x[1]
        # elif near(x[1], Ly):
            # y[0] = x[0]
            # y[1] = x[1] - Ly
        # else:
            # y[0] = -1000
            # y[1] = -1000

# # Sub domain for Periodic boundary condition
# class PeriodicBoundary(SubDomain):

    # # Left boundary is "target domain" G
    # def inside(self, x, on_boundary):
        # xmax = 1.00
        # ymax = 1.00
        # on_left = bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
        # on_bottom = bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS and on_boundary)
        # on_right = bool(x[0] < (xmax + DOLFIN_EPS) and x[0] > (xmax-DOLFIN_EPS) and on_boundary)
        # on_bottom = bool(x[1] < (ymax+DOLFIN_EPS) and x[1] > (ymax-DOLFIN_EPS) and on_boundary)
        # return (on_left or on_right or on_bottom or on_top)

    # def map(self, h, g):
        # xmax = 1.00
        # ymax = 1.00
        # on_left = float(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
        # on_bottom = float(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS and on_boundary)
        # on_right = float(x[0] < (xmax + DOLFIN_EPS) and x[0] > (xmax-DOLFIN_EPS) and on_boundary)
        # on_bottom = float(x[1] < (ymax+DOLFIN_EPS) and x[1] > (ymax-DOLFIN_EPS) and on_boundary)
        # g[0] = h[0] - 1.0
        # g[1] = h[1] - 1.0
        # # if on_g:
            # # g[0] = h[0]
            # # g[1] = h[1] - 1.0
        # # else:
            # # g[0] = h[0]
            # # g[1] = h[1]
            # # print("PBC not applied...")

# Create periodic boundary condition
# pbc = PeriodicDomain()

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

_u_1, _u_2 = u_n.split()
vtkfile_u_1 << (_u_1, 0.0)
vtkfile_u_2 << (_u_2, 0.0)
# plot(u_n2)
# interactive()
sys.exit()

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
  + ((u_2 - u_n2) / k)*v_2*dx + delta*dot(grad(u_2), grad(v_2))*dx + (beta*u_2*(1.0+alpha*r1*u_1*u_2/beta) + u_1*(gamma+r2*u_2))*v_2*dx \
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

