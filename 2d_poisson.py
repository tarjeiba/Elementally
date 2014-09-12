#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
from Integrators import gaussian as gauss
from Meshers import mesh_kit_2D as mesh_kit
import matplotlib.pyplot as plt

###########################
# SET UP PARAMETERS:
##########################

# Poisson function:
def f(x):
    return -1.

# Dirichlet function:
def g(x):
    return 0.

# Neumann function:
def h(x):
    return 1.

##########################
# GENERATE MESH DOMAIN:
##########################

mesh = mesh_kit.quarter_annulus_2D( 0.01, 0.0, np.pi/2., (0.,0.), 1.0, 2.0)

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

#Dirichlet part: degrees of freedom are removed and values of U set.
#Neumann part: Additional contributions to the loading vector.

###########################
#   SOLVE THE SYSTEM:
###########################


