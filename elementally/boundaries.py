import numpy as np






###############################################
##
##          DIRICHLET BOUNDARY:
##
###############################################

def impose_dirichlet(boundary_dict, mesh, 
                     loading_vec, *matrices, time=None):
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

            if time is not None:
                loading_vec[facet[0]] = f(point1,time)
                loading_vec[facet[1]] = f(point2,time)
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


    return loading_vec, res_matrices









