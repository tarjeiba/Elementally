import numpy as np
import numpy.linalg as la

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl


def local_mass_1d():
    pass

def local_mass_2d(local_points):
    """Takes in a meshpy.MeshInfo object and returns a corresponding mass
    matrix.

    INPUT
    local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
    OUTPUT
    M_loc (3,3) np.array
    """

    M_loc = np.zeros((len(local_points), len(local_points)))

    coordinates = np.ones((3, 3))
    coordinates[:, 1:] = local_points

    # Inverting the matrix of coordinates results in the coefficients for the three test
    # functions that are non-zero on this element
    coefficients = la.inv(coordinates)

    for alpha in xrange(3):
        for beta in xrange(3):
            integrand = lambda x: np.inner( coefficients[:,alpha], np.append([1],x))*\
                                  np.inner( coefficients[:,beta], np.append([1],x))
            M_loc[alpha, beta] += gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                                        local_points[2], 4, integrand)
    return M_loc




def local_mass_3d():
    pass
