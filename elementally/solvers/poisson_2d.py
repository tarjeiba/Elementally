#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import time

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl

import assemblers
import meshers
import boundaries

##########################
# BOUNDARY CONDITIONS
##########################
##########################
# SET UP PARAMETERS:
##########################

# Poisson function:
def f(x):
    return 1.

# Dirichlet function:
def g(x):
    return 0.

# Neumann function:
def h(x):
    return 0.

# Boundary dictionary:
boundary_dict = {'dir': {2: g, 3: g}, 'neu': {1: h}}



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

assembly = assemblers.Poisson_2d(mesh, f)
A = assembly.A
b = assembly.b
points = assembly.points
print "Mesh data:"
print " Number of points: ", len(points)
print " Number of elements: ", len(mesh.elements)
print "-----------------------"
print " "

##########################
# BOUNDARY CONDITIONS:
##########################
time1 = time.time()
# Imposing Neumann:
b = boundaries.impose_neumann(boundary_dict, mesh, b,
            nodes, weights)

# Imposing Dirichlet:
A,b = boundaries.impose_dirichlet(boundary_dict, mesh, b, A)

time2 = time.time()
print "Imposing BC: ", time2-time1
###########################
#   SOLVE THE SYSTEM:
###########################
time1 = time.time()
u = la.solve(A, b)
print "Solving linear system: ", time.time() - time1
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

ax.plot_trisurf(points[:,0], points[:,1], u,
                triangles = mesh.elements,
                cmap=cm.jet, linewidth=0.2)

plt.show(1)
