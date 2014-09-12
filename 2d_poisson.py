#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
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

mesh = mesh_kit.quarter_annulus_2D( 0.1 )

points = np.array(mesh.points)
print np.array(mesh.facets)
for face in np.array(mesh.facets):
    print face, mesh.points[face[0]], mesh.points[face[1]]

plt.triplot(points[:,0], points[:,1],np.array(mesh.elements))
plt.show()
##########################
# ASSEMBLY:
##########################

#Needs: 
#   Our understanding of the data structure.

##########################
# BOUNDARY CONDITIONS:
##########################

#Dirichlet part: degrees of freedom are removed and values of U set.
#Neumann part: Additional contributions to the loading vector.

###########################
#   SOLVE THE SYSTEM:
###########################


