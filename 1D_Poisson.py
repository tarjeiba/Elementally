#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrators import gaussian as gauss
from Integrators import GL_quad_functions as gl
import matplotlib as matplot
import matplotlib.pyplot as plt
import scipy
import scipy.linalg

#Set order for quadrature and get nodes and weights:
n = 4
z, w = gl.GL_nodes_and_weights ( n )

# Initiating characteristics for a uniform 1D mesh
x_0 = 0
x_N = 1
num_elements = 100
nodes = np.linspace(x_0, x_N, num_elements + 1)
h = nodes[1] - nodes[0]
f = lambda x: x


# Setting up the boundary conditions
# Dirichlet on the left end point:
u_0 = 1
# Neumann on the right end point:
du_end = -1

# Establishing one possible exact solution
B = du_end + 0.5*x_N ** 2
C = 1.0/6.0 * x_0 ** 3 - B * x_0 + u_0
u_exact = -1./6 * nodes ** 3 + B*nodes + C

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
                                   z,
                                   w)
    b[i+1] += gauss.gaussian_quad_1d(nodes[i], nodes[i+1],
                                   lambda x: ( x - nodes[i] ) / h * f(x),
                                   z,
                                   w)

# MODYFYING MATRIX AND LOADING VECTOR FOR DIRICHLET CONDITION:
A[0,:] = 0
A[0,0] = 1
b[0] = u_0 
# ADDING CONTRIBUTION FROM NEUMANN CONDITION TO LOADING VECTOR:
b[-1] += du_end
# SOLVING RESULTING SYSTEM:
u = scipy.linalg.solve(A,b)

# PLOTTING SOLUTION
# Using plot objects in stead of single commands
# This makes two objects. One image object, which is the entire window that opens,
# and one axes objects, containing the subplots.
# This supposedly greatly simplifies working with the plots further one.
fig, axes = plt.subplots(1,1) # To set a single graph in the figure

plotfontsize = 14
legendfontsize = 12

# Setting the axes object (which is a list if there are multiple plots in one window
# to plot as above

elemental_plot = axes.plot(nodes, u, label="Elemental solution", color="b")
exat_plot = axes.plot(nodes, u_exact, label="An exact solution", color="r")

# Handling the legend of the plot
handles, labels = axes.get_legend_handles_labels()
# inverting their order
axes.legend(handles[::-1], labels[::-1], prop={'size':legendfontsize})

axes.set_xlim(nodes[0], nodes[-1])
axes.set_xlabel("Nodal position (-)", fontsize=plotfontsize, fontweight="normal")
axes.set_ylabel("Value of $u$ (-)", fontsize=plotfontsize, fontweight="normal")

fig.savefig("1d_poisson_solution.png")
plt.show()
