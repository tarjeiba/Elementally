import numpy as np
import numpy.linalg as la

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl


def local_stiffness_1d():
    pass

def local_stiffness_2d(local_points):
    """Takes in a meshpy.MeshInfo object and returns a corresponding stiffness
    matrix.

    INPUT
    local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
    OUTPUT
    A_loc (3,3) np.array
    """

    A_loc = np.zeros((len(local_points), len(local_points)))

    coordinates = np.ones((3, 3))
    coordinates[:, 1:] = local_points

    # Inverting the matrix of coordinates results in the coefficients for the three test
    # functions that are non-zero on this element
    coefficients = la.inv(coordinates)
    area = gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                  local_points[2], 1, lambda x: 1.)

    for alpha in xrange(3):
        for beta in xrange(3):
            A_loc[alpha, beta] += np.inner(coefficients[1:, alpha],
                                           coefficients[1:, beta]) * area
    return A_loc




def local_stiffness_3d():
    pass
