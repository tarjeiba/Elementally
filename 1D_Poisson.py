#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrators import gaussian as gauss
import matplotlib.pyplot as plt
import scipy
import scipy.linalg as linalg

# Initiating characteristics for a uniform 1D mesh
x_0 = 0
x_N = 10
num_elements = 100
nodes = np.linspace(x_0, x_N, num_elements + 1)
h = nodes[1] - nodes[0]
f = lambda x: x
print nodes
# Setting up the boundary conditions
# Dirichlet on the left end point:
u_0 = 0.
# Neumann on the right end point:
du_end = -10.

# ASSEMBLY OF STIFFNESS MATRIX:
# Initialise A-matrix:
A = np.zeros((len(nodes), len(nodes)))
for i in range(num_elements):
    A[i,i] += 1./h
    A[i,i+1] += -1./h
    A[i+1,i] += -1./h
    A[i+1,i+1] += 1./h

# ASSEMBLY OF LOADING VECTOR

# b = np.transpose(np.zeros(len(nodes)))
b = np.zeros(len(nodes))

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
print A
# ADDING CONTRIBUTION FROM NEUMANN CONDITION TO LOADING VECTOR:
b[-1] += du_end
print b
# SOLVING RESULTING SYSTEM:
u = linalg.solve(A, b)

print u


# Finding the exact solution
x_exact = scipy.linspace(x_0, x_N, 100)
u_exact = -1./6 * x_exact ** 3 + ( du_end + 1./2 * x_N ** 2 ) * x_exact + u_0 

# PLOTTING SOLUTION

plt.close("all")

fig = plt.figure()
plt.plot(nodes, u, 'b')
plt.plot(x_exact, u_exact ,'r')
plt.show()

print "Bra"
