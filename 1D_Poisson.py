#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrators import gaussian as gauss
import matplotlib as matplot
import matplotlib.pyplot as pyplot
import scipy

# Initiating characteristics for a uniform 1D mesh
x_0 = 0
x_N = 10
num_elements = 10
nodes = np.linspace(x_0, x_N, num_elements + 1)
h = nodes[1] - nodes[0]
f = lambda x: x
print nodes
# Setting up the boundary conditions
# Dirichlet on the left end point:
u_0 = 1
# Neumann on the right end point:
du_end = -1


# ASSEMBLY OF STIFFNESS MATRIX:
# Initialise A-matrix:
A = np.zeros((len(nodes), len(nodes)))
for i in range(num_elements):
    A[i,i] += 1./h
    A[i,i+1] += -1./h
    A[i+1,i] += -1./h
    A[i+1,i+1] += 1./h

# ASSEMBLY OF LOADING VECTOR

b = np.transpose(np.zeros(len(nodes)))
for i in range(num_elements):
    b[i] += gauss.gaussian_quad_1d(nodes[i], nodes[i+1],
                                   lambda x: ( nodes[i+1] - x ) / h * f(x),
                                   4)
    b[i+1] += gauss.gaussian_quad_1d(nodes[i], nodes[i+1],
                                   lambda x: ( x - nodes[i] ) / h * f(x),
                                   4)
print "Loading vector:"
print b

# MODYFYING MATRIX AND LOADING VECTOR FOR DIRICHLET CONDITION:
A[0,:] = 0
A[0,0] = 1
b[0] = u_0 
# ADDING CONTRIBUTION FROM NEUMANN CONDITION TO LOADING VECTOR:
b[-1] += du_end
print b
# SOLVING RESULTING SYSTEM:
u = scipy.linalg.solve(A,b)

# PLOTTING SOLUTION
pyplot.plot(nodes,u,'b')
u_exact = -1./6 * nodes ** 3 + 4.8*nodes + 1
pyplot.plot(nodes,u_exact,'r')
pyplot.show()
