#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrator import gaussian
import matplotlib as matplot

# Initiating characteristics for a uniform 1D mesh
x_0 = 0
x_N = 10
num_elements = 10
nodes = np.linspace(x_0, x_N, num_elements + 1)
h = nodes[1] - nodes[0]

# ASSEMBLY OF STIFFNESS MATRIX:
# Initialise A-matrix:
A = np.zeros((len(nodes), len(nodes)))
for i in range(num_elements):
    A[i,i] = A[i,i] + 1./h
    A[i,i+1] = A[i,i+1] - 1./h
    A[i+1,i] = A[i+1,i] - 1./h
    A[i+1,i+1] = A[i+1,i+1] + 1./h

print "Nodes: "
print nodes
print "Stiffness matrix: "
print A

# ASSEMBLY OF LOADING VECTOR

b = np.transpose(np.zeros(len(nodes)))

