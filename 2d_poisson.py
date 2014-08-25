#!/usr/bin/env python
#
# This is an initial attempt of an 1D Poisson solver
# written by Trygve and Tarjei, June 26 2014
#

import numpy as np
from Integrators import gaussian as gauss
import matplotlib.pyplot as plt

p1 = np.array([0., 1.])
p2 = np.array([1., 1.])
p3 = np.array([1., 2.])

nq = 4

def f(x):
#    return np.exp(np.sum(x))
    return 1

print gauss.gaussian_quad_2d(p1, p2, p3, nq, f)
