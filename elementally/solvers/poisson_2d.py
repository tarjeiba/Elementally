#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
from meshers import mesh_kit_2d as mesh_kit
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import time

from assemblers import stiffness
from assemblers import loading
from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl

###########################
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

# Order of Gaussian quadrature:
nq = 4
nodes, weights = gl.gl_nodes_and_weights(nq)

##########################
# GENERATE MESH DOMAIN:
##########################

mesh = mesh_kit.quarter_annulus_2d( 0.01, 0.0, 2.*np.pi/2., (0.,0.), 1.0, 2.0)

##########################
# ASSEMBLY:
##########################

points = np.array(mesh.points)
A = np.zeros((len(mesh.points), len(mesh.points)))
b = np.zeros(len(mesh.points))

time1 = time.time()
for element in mesh.elements:
    A[np.ix_(element, element)] += stiffness.local_stiffness_2d(points[element])
    b[element] += loading.local_loading_2d(points[element], f)

time2 = time.time()
print "Mesh data:"
print " Number of points: ", len(points)
print " Number of elements: ", len(mesh.elements)
print "-----------------------"
print " "
print "Assembly time: ", time2-time1


##########################
# BOUNDARY CONDITIONS:
##########################
time1 = time.time()
updated_facets = mesh_kit.facet_orientation(mesh.facets, mesh.elements, points)
# Neumann boundary:
for i, facet in enumerate(updated_facets):
    if (mesh.facet_markers[i] == 1):
        p1 = points[facet[0],:]
        p2 = points[facet[1],:]
        coeffs = la.inv( np.vstack( (p1,p2) ) )
        #Adding the contributions:
        b[facet[0]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,0],x) )
        b[facet[1]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,1],x) )

for i, facet in enumerate(updated_facets):
    if (mesh.facet_markers[i] == 2):
        A[facet[0],:] = 0
        A[facet[1],:] = 0
        A[facet[0], facet[0]] = 1
        A[facet[1],facet[1]] = 1
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

plt.show(1)

#print mesh.facets
#for i, facet in enumerate(mesh.facets):
#    print facet, mesh.facet_markers[i]
