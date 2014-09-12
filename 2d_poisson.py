#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrators import gaussian as gauss
import matplotlib.pyplot as plt

###########################
#SET UP PARAMETERS:
##########################


##########################
#GENERATE MESH DOMAIN:
##########################

##########################
# ASSEMBLY:
##########################

#Needs: Line integral (for loading vector).
#   2D quadrature for stiffness matrix and load vector.
#   Our understanding of the data structure.

##########################
# BOUNDARY CONDITIONS:
##########################

#Dirichlet part: degrees of freedom are removed and values of U set.
#Neumann part: Additional contributions to the loading vector.

###########################
#   SOLVE THE SYSTEM:
###########################


