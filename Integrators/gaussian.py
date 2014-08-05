# This is the module file for the gaussian integrator
# Created June 26 2014

import numpy as np

def gaussian_quad_1d(a, b, f, n):
    """One-dimensional integrator using Gaussian quadrature.
    INPUT:
        a,b: End points of interval integration is taken over.
        f: The integrand passed to the function as a function handle (lambda).
        n: Number of quadrature nodes. Per June 2014 n need to be in the set {1,2,3,4}
    OUTPUT:
        Returns the approximated integral as a real.
    Typical usage:
        Example: Suppose you want to integrate exp(x) over the interval (1,2) using 3 nodes.
        Simply call
            gaussian_quad_1d(1, 2, lambda x: np.exp(x), 3)
    """
    
    #Getting out nodes and weights for the referance interval (-1,1):
    z,w = nodes_and_weights_1d(n)
    #Scaling nodes and weights for relevant interval (a,b):
    z = (b + a)/2. + (b - a)/2. * z
    w = (b - a)/2. * w

    #Returning approximate integral:    
    return np.inner(w, f(z))
    
    pass

def gaussian_quad_2d(p1, p2, p3, nq, f):
    """Two-dimensional integrator using Gaussion quadrature.
    INPUT:
        p1, p2, p3: Corner points of triangle.
        nq: Number of quadrature points.
        f: The integrator passed as a function handle (lambda).
    OUTPUT:
        Returns the approximated integral as a real.
    Typical usage:
    HOLD
    """
    nodes, weights =  nodes_and_weights_2d(nq)
    
    coordinates = np.vstack ( (p1, p2, p3) )
    coordinates = np.dot ( nodes, coordinates )

    
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

