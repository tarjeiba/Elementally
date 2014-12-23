#!/usr/bin/env python
#
# This is an initial attempt of an 2D Poisson solver
# written by Trygve and Tarjei, September 2014
#

import numpy as np
import numpy.linalg as la
from meshers import mesh_kit_2d as mesh_kit
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

from assemblers import stiffness
from assemblers import loading
from assemblers import mass
from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl

###########################
# SET UP PARAMETERS:
##########################

# Time-step
dt = 1.0E-4

# Start time
t_start = 0.0

# stop time
t_stop = 1.0

# Directioness of Euler
theta = 0.5

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
M = np.zeros((len(mesh.points), len(mesh.points)))
A = np.zeros((len(mesh.points), len(mesh.points)))
b = np.zeros(len(mesh.points))

for element in mesh.elements:
    M[np.ix_(element, element)] += stiffness.local_stiffness_2d(points[element])
    A[np.ix_(element, element)] += stiffness.local_stiffness_2d(points[element])
    b[element] += loading.local_loading_2d(points[element], f)

W1 = M + dt * (1 - theta) * A
W2 = M - dt * theta * A
    
##########################
# INITIAL CONDITIONS
##########################

c = np.zeros(len(mesh.points))


##########################
# BOUNDARY CONDITIONS:
##########################

updated_facets = mesh_kit.facet_orientation(mesh.facets, mesh.elements, points)
# Neumann boundary:
for i, facet in enumerate(updated_facets):
    if (mesh.facet_markers[i] == 1):
        p1 = points[facet[0],:]
        p2 = points[facet[1],:]
        coeffs = la.inv(np.vstack((p1,p2)))
        #Adding the contributions:
        b[facet[0]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,0],x))
        b[facet[1]] += gauss.gaussian_line(p1, p2, nodes, weights,
                lambda x: h(x) * np.inner(coeffs[:,1],x))


# Dirichlet boundary:
for i, facet in enumerate(updated_facets):
    if (mesh.facet_markers[i] == 2):
        W1[facet[0],:] = 0
        W1[facet[1],:] = 0
        W2[facet[0],:] = 0
        W2[facet[1],:] = 0
        W1[facet[0], facet[0]] = 1
        W1[facet[1],facet[1]] = 1
        b[facet[0]] = g(points[facet[0],:])
        b[facet[1]] = g(points[facet[1],:])

#########################################################
#       SOLVING AND ANIMATING:
#########################################################
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = None
t = t_start
while (t<t_stop):
    oldcol = surf

    c = la.solve(W1, np.dot(W2, c) + b)
    surf = ax.plot_trisurf(points[:,0], points[:,1], c,
                            triangles = mesh.elements,
                            cmap=cm.jet, linewidth=0.2)


    if oldcol is not None:
        ax.collections.remove(oldcol)

    plt.pause(0.01)
    t += dt



def plot_init():
    pass

def plot_animate(ax):
    pass

def step(t, dt, c, W1, W2):
    c_0 = np.copy(c)
    c = la.solve(W1, np.dot(W2, c_0) + b)
    return c










plt.show(1)

#print mesh.facets
#for i, facet in enumerate(mesh.facets):
#    print facet, mesh.facet_markers[i]
