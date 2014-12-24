import numpy as np
import scipy.linalg as la
from integrators import gaussian as gauss






###############################################
##
##          DIRICHLET BOUNDARY:
##
###############################################

def impose_dirichlet(boundary_dict, mesh, 
                     loading_vec, *matrices,**kwargs):
    """Function for imposing the boundary Dirichlet boundary,
    conditions.
    INPUT:
        boundary_dict: Dictionary containing facet_markers and
            corresponding functions.
        mesh: ElementallyMeshInfo, here we need facets, facet_markers,
            and points.
        loading_vec: Loading vector for which the Dirichlet is imposed.
        *matrices: Unknown number of matrices for which the rows are
            to be zeroed out. For the first of these the diagonal
            element is set to unity.
        time=None: Allows for time dependent functions if time is set
            to something different from None.
    OUTPUT:
        Updated matrices and loading vector.
    """
    if 'time' in kwargs.keys():
        t = kwargs['time']
    else:
        t = None

    res_matrices = []
    for matrix in matrices:
        res_matrices.append(matrix)
    # Let's start by iterating over each facet.
    for i, facet in enumerate(mesh.facets):
        # Check if current facet is indeed a Dirichlet facet:
        if mesh.facet_markers[i] in boundary_dict['dir'].keys():
            # Get points on facet:
            point1 = np.array(mesh.points[facet[0]])
            point2 = np.array(mesh.points[facet[1]])

            # Fix element in loading vector:
            f = boundary_dict['dir'][mesh.facet_markers[i]]

            if t is not None:
                loading_vec[facet[0]] = f(point1,t)
                loading_vec[facet[1]] = f(point2,t)
            else:
                loading_vec[facet[0]] = f(point1)
                loading_vec[facet[1]] = f(point2)

            #Go through all matrices:
            for j, matrix in enumerate(res_matrices):
                # Zero out row:
                matrix[facet[0]] = 0.
                matrix[facet[1]] = 0.

                # If first row; unity diagonal element:
                if (j==0):
                    matrix[facet[0],facet[0]] = 1.
                    matrix[facet[1],facet[1]] = 1.


    if (len(res_matrices)==1):
        return res_matrices[0], loading_vec
    else:
        return res_matrices, loading_vec


################################################################
##
##                  NEUMANN BOUNDARY
##
################################################################

def impose_neumann(boundary_dict,mesh, loading_vec,
                   nodes, weights, time=None):
    """Function for imposing Neumann boundary conditions.
    INPUT:
        boundary_dict: Dictionary with boundary data.
        mesh: ElementallyMeshInfo.
        loading_vec: Loading vector.
        nodes, weights: Used for numerical quadrature.
        time(defalt=None): If time dependent problem.
    OUTPUT:
        Updated loading vector.
    """
    # We start by iterating over facets:
    for i, facet in enumerate(mesh.facets):
        # Check if the facet is Neumann:
        if mesh.facet_markers[i] in boundary_dict['neu'].keys():
            #Get points:
            p1 = np.array(mesh.points[facet[0]])
            p2 = np.array(mesh.points[facet[1]])
            coeffs = la.inv( np.vstack( (p1,p2) ) )
            # Create proto-integrand:
            if time is not None:
                f = lambda x: \
                    boundary_dict['neu'][mesh.facet_markers[i]](x,time)
            else:
                f = lambda x: \
                    boundary_dict['neu'][mesh.facet_markers[i]](x)

            # Add the contributions:
            loading_vec[facet[0]] += gauss.gaussian_line(p1,p2,
                nodes, weights,
                lambda x: f(x) * np.inner(coeffs[:,0],x) )

            loading_vec[facet[1]] += gauss.gaussian_line(p1,p2,
                nodes, weights,
                lambda x: f(x) * np.inner(coeffs[:,1],x) )

    return loading_vec




