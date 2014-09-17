# This is the module file for the gaussian integrator
# Created June 26 2014

import numpy as np

def gaussian_quad_1d(a, b, f, z, w):
    """One-dimensional integrator using Gaussian quadrature.
    INPUT:
        a,b: End points of interval integration is taken over.
        f: The integrand passed to the function as a function handle (lambda).
        z: Quadrature nodes.
        w: Quadrature weights.
    OUTPUT:
        Returns the approximated integral as a real.
    Typical usage:
        Example: Suppose you want to integrate exp(x) over the interval (1,2) using 3 nodes.
        Simply call
            gaussian_quad_1d(1, 2, lambda x: np.exp(x), 3)
    """
    
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

    weights = np.linalg.det( np.vstack( (p2-p1,p3-p1) ) ) * weights
    return 0.5*np.inner(weights, np.apply_along_axis(f, 1, coordinates) )

    
    pass

def nodes_and_weights_1d(nq):
    """Takes the number of nodes (from 1 to 4) and returns
    the nodal positions and their respective weights as
    two distinct np.arrays.
    
    Per June 28 2014, this has been set to be improved.
    """
    z = np.array([])
    w = np.array([])
    if nq == 1:
        z = np.append(z, 0)
        w = np.append(w, 2)
        return z, w
    elif nq == 2:
        z = np.append(z, [-1 * np.sqrt(1./3.), np.sqrt(1./3.)])
        w = np.append(w, [1., 1.])
        return z, w
    elif nq == 3:
        z = np.append(z, [-1 * np.sqrt(3./5), 0, np.sqrt(3./5)])
        w = np.append(w, [5./9., 8./9., 5./9.])
        return z, w
    elif nq == 4:
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
        return 1

def nodes_and_weights_2d(nq):
    """Takes the number of quadrature points (1, 2, or 3) and returns
    the nodal position and their respective weights as two distinct
    np.arrays.

    This is to be improved to accept an arbitrary number of quadrature points.
    """

    if nq == 1:
        z = np.array([[1./3, 1./3, 1./3]])
        w = np.array([1.])
        return z, w

    elif nq == 3:
        z = np.array([[1./2, 1./2, 0],
                      [1./2, 0, 1./2],
                      [0, 1./2, 1./2]])
        w = np.array([1./3, 1./3, 1./3])
        return z, w

    elif nq == 4:
        z = np.array([[1./3, 1./3, 1./3],
                      [3./5, 1./5, 1./5],
                      [1./5, 3./5, 1./5],
                      [1./5, 1./5, 3./5]])
        w = np.array([-9./16, 25./48, 25./48, 25./48])
        return z, w

    else:
        print "Please use a valid number of nodes (1, 3, or 4)."
        return 1


def gaussian_line(start, end, nodes, weights, f):
    """Calculates the numerical line integral of the function f from start to end, using pre-defined
    nodes and weights.
    INPUT:
        start, end: one-dimensional vectors in R2
        nodes, weights: one-dimensional vector of arbitrary length
        f: function handle returning a scalar
    OUTPUT:
        Approximate integral as a real
    """
    # Finding the world coordinate evaluation points
    world_nodes = np.outer(nodes, (end - start)/2.) + (end + start) / 2.
    world_weights = np.sqrt(np.inner(end - start, end - start)) * weights / 2.

    # Applying the function f along each row of the wolrd_nodes matrix
    evals = np.apply_along_axis(f, 1, world_nodes)

    return np.inner(world_weights, evals)

