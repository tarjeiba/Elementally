#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl
from meshers import mesh_kit_2d as mesh_kit
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

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

for element in mesh.elements:
    # i, j, and k are the indices of the points defining the triangular element
    coord_i, coord_j, coord_k = points[element, :]
    coordinates = np.ones((3, 3))
    coordinates[:, 1:] = np.vstack((coord_i, coord_j, coord_k))
    # Inverting the matrix of coordinates results in the coefficients for the three test
    # functions that are non-zero on this element
    coefficients = la.inv(coordinates)
    area = gauss.gaussian_quad_2d(coord_i, coord_j, coord_k, 1,lambda x: 1.)
    for alpha in xrange(3):
        i = element[alpha]
        integrand = lambda x: f(x) * np.inner(coefficients[:, alpha], np.append([1], x))
        b[i] += gauss.gaussian_quad_2d(coord_i, coord_j, coord_k, 4, integrand)
        for beta in xrange(3):
            j = element[beta]
            A[i, j] += np.inner(coefficients[1:, alpha], coefficients[1:, beta]) * area

##########################
# BOUNDARY CONDITIONS:
##########################

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

###########################
#   SOLVE THE SYSTEM:
###########################

u = la.solve(A, b)
fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

ax.plot_trisurf(points[:,0], points[:,1], u, triangles = mesh.elements, cmap=cm.jet, linewidth=0.2)

plt.show(1)

#print mesh.facets
#for i, facet in enumerate(mesh.facets):
#    print facet, mesh.facet_markers[i]
