#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
import matplotlib as matplot

# Initiating charactestics for a uniform 1D mesh
x_0 = 0
x_N = 10
num_nodes = 10
nodes = np.linspace(x_0, x_N, num_nodes)

print nodes
