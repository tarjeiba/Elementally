import numpy as np
import numpy.linalg as la

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl


def local_loading_1d():
    pass

def local_loading_2d(local_points, load_func):
    """Takes in a meshpy.MeshInfo object and returns a corresponding loading
    vector.

    INPUT
    local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
    load_func function handle for scalar returning loading function
    OUTPUT
    b_loc loading vector of length 3
    """

    b_loc = np.zeros(len(local_points))

    coordinates = np.ones((3, 3))
    coordinates[:, 1:] = local_points

    # Inverting the matrix of coordinates results in the coefficients for
    # the three test functions that are non-zero on this element
    coefficients = la.inv(coordinates)

    for alpha in xrange(3):
        integrand = lambda x: load_func(x) * np.inner(coefficients[:, alpha],
                                              np.append([1], x))
        b_loc[alpha] += gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                               local_points[2], 4, integrand)
    return b_loc

def stiffness_3d():
    pass
