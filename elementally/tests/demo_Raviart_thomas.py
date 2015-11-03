#!/usr/bin/env python
#
# Demo of creating a Raviart-Thomas function space, and we plot it on
# the reference element.
##----------------------------------##

import numpy as np
import matplotlib.pyplot as plt

import meshers
import function_spaces

#######################################

# Create mesh (Needed to construct function space)
mesh = meshers.unit_square_2d(2,2,\
        generate_facets=True,\
        generate_neighbors=True)

# Create Raviart-Thomas function space:
RTh = function_spaces.RaviartThomas(mesh)

# Create meshgrid to plot on:
N = 10
x = np.linspace(0., 1., N)

X, Y = np.meshgrid(x,x, indexing="ij")

# Plot basis functions:
k = 2
U = np.zeros ( (N, N) )
V = np.zeros ( (N, N) )
for i in range(N):
  for j in range(N):
    x = np.array( (X[i,j], Y[i,j]) )
    U[i,j], V[i,j] = RTh.funcs[k](x)

plt.figure()
plt.quiver(X,Y,U,V)
plt.show()
