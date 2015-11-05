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
N = 40
mesh = meshers.unit_square_2d(N, N,\
          generate_faces=True,\
          generate_normals=True)

print "Finished creating mesh, starting on edge connectivity"
print "Total number of edges: ", len(mesh.faces)
t1 = time.time()
# Need edge connectivity:
mesh.set_edges_of_elements()

t2 = time.time()
print "Time: ", t2-t1

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

plt.spy(A)
plt.show()


########################
##    POST-PROCESSING:
########################
points = np.array(mesh.points)
elements = np.array(mesh.elements)

fig = plt.figure(1)
ax = fig.gca()

ax.triplot(points[:,0], points[:,1],\
          triangles=elements)

plt.show()
