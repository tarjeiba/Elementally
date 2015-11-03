#!/usr/bin/env python
#
# Demo of creating a Raviart-Thomas function space, and we plot it on
# the reference element.
##----------------------------------##

import numpy as np
import matplotlib.pyplot as plt

import meshers
import problem_class

#######################################

# Create mesh (Needed to construct function space)
mesh = meshers.unit_square_2d(2,2,\
        generate_facets=True,\
        generate_neighbors=True)

# Create Raviart-Thomas function space:
RTh = RaviartThomas(mesh)
