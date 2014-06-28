# This is the module file for the gaussian integrator
# Created June 26 2014

import numpy as np

def gaussian_quad_1d():
    """
    One-dimensional integrator using Gaussian quadrature.
    Typical usage:


    """

    pass

def nodes_and_weights_1d(n):
    """
    Takes the number of nodes (from 1 to 4) and returns
    the nodal positions and their respective weights as
    two distinct np.arrays.
    
    Per June 28 2014, this has been set to be improved.
    """
    z = np.array([])
    w = np.array([])
    if n == 1:
        z = np.append(z, 0)
        w = np.append(w, 2)
        return z, w
    elif n == 2:
        z = np.append(z, [-1 * np.sqrt(1./3.), np.sqrt(1./3.)])
        w = np.append(w, [1., 1.])
        return z, w
    elif n == 3:
        z = np.append(z, [-1 * np.sqrt(3./5), 0, np.sqrt(3./5)])
        w = np.append(w, [5./9., 8./9., 5./9.])
        return z, w
    elif n == 4:
        z = np.append(z, [-1 * np.sqrt((3. + 2. * np.sqrt(6./5.)) / 7.),
                          -1 * np.sqrt((3. - 2. * np.sqrt(6./5.)) / 7.),
                          np.sqrt((3. - 2. * np.sqrt(6./5.)) / 7.),
                          np.sqrt((3. + 2. * np.sqrt(6./5.)) / 7.)])
        w = np.append(w, [(18. - np.sqrt(30.)) / 36.,
                         (18. + np.sqrt(30.)) / 36.,
                         (18. + np.sqrt(30.)) / 36.,
                         (18. - np.sqrt(30.)) / 36.])
        return z, w
    else:
        print "Please use a valid number of nodes"
        return 0

