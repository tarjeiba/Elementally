#!/usr/bin/env python
#
# An example of using Mixed Finite elements for
# solving Poisson equation, November 2015.
#

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl

import assemblers
import meshers
import boundaries

import time
####################
##  PARAMETERS:
####################

# Loading function:
def f(x):
    return 1.

# Dirichlet:
def g(x):
    return 0.

# Neumann:
def h(x):
    return 0.

# Order of gaussian quadrature:
nq = 4
nodes, weights = gl.gl_nodes_and_weights(nq)

# Boundary dictionary:
boundary_dict = {'dir': {1: g, 3: g}, 'neu': {2: h}}

########################
##    GENERATE MESH
########################
t1 = time.time()
N = 20
mesh = meshers.unit_square_2d(N, N,\
          generate_faces=True,\
          generate_normals=True)

print "Finished creating mesh, starting on edge connectivity"
print "Total number of edges: ", len(mesh.faces)
print "Total number of triangles: ", len(mesh.elements)
# Need edge connectivity:
mesh.set_edges_of_elements()

t2 = time.time()
print "Meshing time: ", t2-t1

########################
##    ASSEMBLY:
########################
assembly = assemblers.MixedPoisson_2d(mesh, f)
A = assembly.A
b = assembly.b
B = assembly.B

# Check sparsity pattern:
A = sp.csr_matrix(A)
B = sp.csr_matrix(B)

#Check assembly time:
t3 = time.time()
print "Assembly time: ", t3-t2

# Assemble total matrix:
C = sp.bmat( [[A, B.T],\
              [B, None]],\
              format='lil')

##########################
##  BOUNDARY CONDITIONS:
##########################
t1 = time.time()
# Impose Dirichlet weakly:
G = np.zeros(assembly.dofs_sig)
boundaries.impose_dirichlet_mixed_poisson(boundary_dict,\
            mesh, G,\
            nodes, weights)


# Total right hand side:
F = np.concatenate( (G,b) )

# Impose Neumann strongly:
boundaries.impose_neumann_mixed_poisson(boundary_dict, mesh,\
                    nodes, weights,\
                    F, C)
# Set to CSR-format:
C.tocsr()
t2 = time.time()
print "Time to impose BC: ", t2-t1
#########################
##    SOLVE:
#########################
# We try to use a iterative solver:

t1 = time.time()
# NOTE: spla.minres, although fast, is as of yet
# not properly tested by SciPy.
# Seems like lgmres is a viable option concerning
# speed for now. Another advantage here is that symmetry
# is not assumed.
X, info = spla.lgmres(C,F)
t2 = time.time()
print "Time to solve: ", t2-t1
if not info:
    print "Iterative solver exited successfully."
else:
    print "Iterative solver had an issue. Exit status: ", info

# Split solutions:
sigma = X[:assembly.dofs_sig]
u = X[assembly.dofs_sig:]

########################
##    POST-PROCESSING:
########################
points = np.array(mesh.points)
elements = np.array(mesh.elements)

fig = plt.figure(1)
ax = fig.gca()

# Plot u:
ax.tripcolor(points[:,0], points[:,1], facecolors=u,\
          triangles=elements,\
          edgecolors='k')
plt.show()
