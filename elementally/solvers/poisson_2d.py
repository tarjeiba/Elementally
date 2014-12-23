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

##########################
# BOUNDARY CONDITIONS
##########################

boundaries = {
    'dir': [(1, g1),
            (2, g2)],
    'neu': [(3, h1)]
    }


##########################
# SET UP PARAMETERS:
##########################

# Poisson function:
def f(x):
    return 1.

# Dirichlet function:
def g(x):
    return 0.

def g1(x):
    return 0.

def g2(x, t):
    return 1. * np.sin(4 * np.pi * t / t_stop)
# Neumann function:
def h(x):
    return 0.

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
# Neumann boundary:
for i, facet in enumerate(mesh.facets):
    if (mesh.facet_markers[i] == 1):
        p1 = points[facet[0],:]
        p2 = points[facet[1],:]
        coeffs = la.inv( np.vstack( (p1,p2) ) )
        #Adding the contributions:
        b[facet[0]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,0],x) )
        b[facet[1]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,1],x) )

for i, facet in enumerate(mesh.facets):
    if (mesh.facet_markers[i] == 2):
        A[facet[0],:] = 0
        A[facet[1],:] = 0
        A[facet[0], facet[0]] = 1
        A[facet[1], facet[1]] = 1
        b[facet[0]] = g(points[facet[0],:])
        b[facet[1]] = g(points[facet[1],:])

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
