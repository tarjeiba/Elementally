#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl

import assemblers
import meshers
import boundaries
###########################
# SET UP PARAMETERS:
##########################

# Time-step
dt = 1.0E-1

# Start time
t_start = 0.0

# stop time
t_stop = 20

# Directioness of Euler
theta = 0.5
##############################
# Poisson function:
def f(x):
    return 0.

# Dirichlet functions:
def g1(x,t):
    return 0.

def g2(x, t):
    return 1. * np.sin(8 * np.pi * t / t_stop)

# Neumann function:
def h(x,t):
    return 0.

boundary_dict = {'dir': {3: g2}, 'neu': {1: h, 2: g1}}


# Order of Gaussian quadrature:
nq = 4
nodes, weights = gl.gl_nodes_and_weights(nq)

##########################
# GENERATE MESH DOMAIN:
##########################

mesh = meshers.quarter_annulus_2d( 0.01, 0.0, 2.*np.pi/2., (0.,0.), 1.0, 2.0)

##########################
# ASSEMBLY:
##########################

assembly = assemblers.Heat_Equation_2d(mesh,f)
A = assembly.A
M = assembly.M
# Constituents of the loading vector
loading_base = assembly.b
loading_bc = np.zeros(len(mesh.points))
loading_bc_old = np.copy(loading_bc)
loading_dir = np.copy(loading_bc)

points = assembly.points


W1 = M + dt * (1 - theta) * A
W2 = M - dt * theta * A
    
##########################
# INITIAL CONDITIONS
##########################

c = np.zeros(len(mesh.points))

##########################
# BOUNDARY CONDITIONS:
##########################

# Neumann boundary:
loading_bc = boundaries.impose_neumann(boundary_dict, mesh, loading_bc,
                nodes, weights, time=t_start)


# Dirichlet boundary:
[W1, W2], loading_bc = boundaries.impose_dirichlet(boundary_dict, mesh,
                loading_bc, W1, W2, time=t_start)

# Also need to zero out some rows of the loading_base vector.
loading_base = boundaries.impose_dirichlet(boundary_dict, mesh,
                loading_base, time=t_start)


#########################################################
#       SOLVING AND ANIMATING:
#########################################################

fig = plt.figure()
# Set colormap:
colormap = cm.get_cmap('bone')


ax = fig.gca(projection='3d')

ax.set_xlim3d(left=-2., right=2.)
ax.set_ylim3d(bottom=0., top=2.)
ax.set_zlim3d(bottom=-1., top=1.)

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$c$')
ax.set_title('Numerical solution')

surf = ax.plot_trisurf(points[:,0], points[:,1], c,
                       triangles = mesh.elements,
                       cmap=colormap, vmin=-1.,
                       vmax= 1., linewidth=0.2)
t = t_start
while (t<t_stop):
    # Update time:
    t += dt
    # Update old surface:
    oldcol = surf
    # Update old BC loading:
    loading_bc_old = loading_bc

    # Construct new loading_bc:
    loading_bc = np.zeros(len(mesh.points))
    loading_bc = boundaries.impose_neumann(boundary_dict, mesh,
                    loading_bc, nodes, weights, time=t)


    loading_dir = boundaries.impose_dirichlet(boundary_dict, mesh,
                    loading_dir, time=t)
    # create time-loading (just a temp):
    temp = loading_base + (1-theta)*loading_bc + theta*loading_bc_old

    c = la.solve(W1, np.dot(W2, c) + dt*temp + loading_dir)
    surf = ax.plot_trisurf(points[:,0], points[:,1], c,
                            triangles = mesh.elements,
                            cmap=colormap, vmin=-1.,
                            vmax=1., linewidth=0.2)


    if oldcol is not None:
        ax.collections.remove(oldcol)

    plt.pause(0.0001)

plt.show(1)

