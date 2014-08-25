#This file contains functions used for getting nodes and weights
#to be used in Gaussian quadratures.
import numpy as np

def eval_Legendre ( x, n ):
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
    P[1] = x

    for i in range(1,n-1):
        P[i+1] = ( (2*i+1) * x * P[i] - i * P[i-1] ) / float(i+1)

    return P[n]

def eval_Legendre_diff ( x, n ):
    """
