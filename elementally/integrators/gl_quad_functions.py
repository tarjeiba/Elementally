#This file contains functions used for getting nodes and weights
#to be used in Gaussian quadratures.
import numpy as np
import scipy.linalg as la

def eval_legendre ( x, n ):
    """This function evaluates the Legendre polynomial
    of order n at the point x.
    Uses the recursion relation for evaluating.
    INPUT:
        x: Real - point where Legendre polynomial is to be evaluated.
        n: integer - order of polynomial to be evaluated.
    OUTPUT:
        P_n(x): Real.
    """
    P = np.zeros(n+1)
    P[0] = 1.0
    if (n>=1):
        P[1] = x
        for i in range(1,n):
            P[i+1] = ( (2*i+1) * x * P[i] - i * P[i-1] ) / float(i+1)
    
    return P[n]

def eval_legendre_diff ( x, n ):
    """Function for evaluating the derivative of the Legendre
    polynomial of order n.
    INPUT:
        x: Real - point to be evaluated.
        n: Integer - order of Legendre polynomial.
    OUTPUT:
        P_n'(x): Real.
    """
    return n*x/(x ** 2 -1 ) * eval_legendre( x, n ) -n/(x ** 2 - 1) * eval_legendre( x, n-1 )

def gl_nodes_and_weights ( n ):
    """Function for finding nodes and weights
    for n-point Gaussian quadrature.
    INPUT:
        n: Integer - how many points to be used in the ensuing quadratures.
    OUTPUT:
        z: array of points.
        w: array of weights.
    """

    z = np.zeros( n )
    w = np.zeros( n )
    if (n == 1):
        z[0] = 0.0
    else:
        A = np.zeros( (n,n) )
        A[0,1] = 1
        for i in range(1,n-1):
            A[i,i-1] = i / float(2*i+1)
            A[i,i+1] = (i+1) / float(2*i+1)

        A[-1,-2] = (n-1) / float(2*n-1)
        z = np.sort( la.eig(A, right=False) ).real

    for i in range(n):
        w[i] = 2. / ( (1-z[i] ** 2) * (eval_legendre_diff( z[i], n )**2) )

    return z, w
