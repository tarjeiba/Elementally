import numpy as np
import numpy.linalg as la
import scipy.sparse as sp


from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl
import meshers

##############################################################
##
##              BASE ASSEMBLER CLASSES:
##
##############################################################
class Assembly(object):

    def __init__(self):
        pass


class Assembly_1d(Assembly):

    def __init__(self, mesh):
        super(Assembly_1d, self).__init__()

        if self.equation is not None:
            print "Assembly attempted without providing equation."
            return 1

    def local_mass(self):
        pass

    def local_stiffness(self):
        pass

    def local_loading(self):
        pass


class Assembly_2d(Assembly):

    def __init__(self, mesh):
        super(Assembly_2d, self).__init__()
        self.points = np.array(mesh.points)

    def local_mass(self, local_points, local_coeffs):
        """Takes in a meshpy.MeshInfo object and returns a corresponding mass
        matrix.

        INPUT
        local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
        local_coeffs: Coefficients for local basis functions.
        OUTPUT
        M_loc (3,3) np.array
        """

        M_loc = np.zeros((len(local_points), len(local_points)))


        # Inverting the matrix of coordinates results in the coefficients for the three test
        # functions that are non-zero on this element

        for alpha in xrange(3):
            for beta in xrange(3):
                integrand = lambda x: np.inner( local_coeffs[:,alpha], np.append([1],x))*\
                                      np.inner( local_coeffs[:,beta], np.append([1],x))
                M_loc[alpha, beta] += gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                                            local_points[2], 4, integrand)
        return M_loc

    def local_stiffness(self, local_points, local_coeffs):
        """Takes in a meshpy.MeshInfo object and returns a corresponding stiffness
        matrix.

        INPUT
        local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
        OUTPUT
        A_loc (3,3) np.array
        """

        A_loc = np.zeros((len(local_points), len(local_points)))


        # Inverting the matrix of coordinates results in the coefficients for the three test
        # functions that are non-zero on this element
        area = gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                      local_points[2], 1, lambda x: 1.)

        for alpha in xrange(3):
            for beta in xrange(3):
                A_loc[alpha, beta] += np.inner(local_coeffs[1:, alpha],
                                               local_coeffs[1:, beta]) * area
        return A_loc

    def local_loading(self, local_points, load_func, local_coeffs):
        """Takes in a meshpy.MeshInfo object and returns a corresponding loading
        vector.

        INPUT
        local_points (3,2) np.array [[x1,y1], [x2,y2], [x3,y3]]
        load_func function handle for scalar returning loading function
        OUTPUT
        b_loc loading vector of length 3
        """

        b_loc = np.zeros(len(local_points))

        for alpha in xrange(3):
            integrand = lambda x: load_func(x) * np.inner(local_coeffs[:, alpha],
                                                  np.append([1], x))
            b_loc[alpha] += gauss.gaussian_quad_2d(local_points[0], local_points[1],
                                                   local_points[2], 4, integrand)
        return b_loc




class Assembly_3d(Assembly):

    def local_mass(self):
        pass

    def local_stiffness(self):
        pass

    def stiffness(self):
        pass

###################################################
##
##            PROBLEM SPECIFIC ASSEMBLERS:
##
###################################################

##    POISSON 2D    ##
class Poisson_2d(Assembly_2d):

    def __init__(self, mesh, load_func):
        super(Poisson_2d, self).__init__(mesh)
        self.A = np.zeros((len(self.points), len(self.points)))
        self.b = np.zeros(len(self.points))
        self.load_func = load_func
        self.assembler(mesh)

    def assembler(self, mesh):
        for element in mesh.elements:
            local_points = np.array(self.points[element])
            coords = np.ones( (3,3) )
            coords[:,1:] = local_points
            coeffs = la.inv(coords)
            self.A[np.ix_(element, element)] += self.local_stiffness(local_points, coeffs)
            self.b[element] += self.local_loading(local_points, self.load_func, coeffs)

##    HEAT EQUATION 2D    ##
class Heat_Equation_2d(Assembly_2d):

    def __init__(self,mesh, load_func):
        super(Heat_Equation_2d, self).__init__(mesh)
        self.A = np.zeros((len(self.points), len(self.points)))
        self.M = np.zeros((len(self.points), len(self.points)))
        self.b = np.zeros(len(self.points))
        self.load_func = load_func
        self.assembler(mesh)

    def assembler(self, mesh):
        for element in mesh.elements:
            # Coefficients for the basis functions are found by
            # inverting the coordinate matrix.
            local_points = np.array(self.points[element])
            coords = np.ones( (3,3) )
            coords[:,1:] = local_points
            coeffs = la.inv(coords)
            
            self.A[np.ix_(element, element)] += self.local_stiffness(local_points, coeffs)
            self.M[np.ix_(element, element)] += self.local_mass(local_points, coeffs)
            self.b[element] += self.local_loading(local_points, self.load_func, coeffs)


##    MIXED POISSON   ##
########################
# Here we come to a very different case than before.
# First, we will now encounter two new types of elements,
# different from the linear elements we've seen so far
# (basis function defined by value at triangle vertices).
# In addition, this is our first encounter with a mixed element method
# (system of PDEs).
# We consider, yet again Poisson equation, but now both u and sig = grad(u)
# are primary unknowns. Different from before is that u will now be sought in
# the space of piecewise constants over each triangle, and sig will be sought in 
# the 1st degree Raviart-Thomas space; piecewise linear vector fields with continuous
# flux across triangle edges.
class MixedPoisson_2d(Assembly_2d):

    # Initializer:
    def __init__(self, mesh, load_func):
        # Input mesh needs to have all faces with face markers defined
        # and properly set.
        super(MixedPoisson_2d, self).__init__(mesh)
        # Number of degrees of freedom:
        self.points = np.array(mesh.points)
        self.dofs_u = len(mesh.elements)
        self.dofs_sig = len(mesh.faces)
        
        # We use sparse matrices, and the lil_matrix structure, since the
        # sparsity structure is going to change as we iterate over elements.
        self.A = sp.lil_matrix( (self.dofs_sig, self.dofs_sig) )
        self.B = sp.lil_matrix( (self.dofs_u, self.dofs_sig) )
        self.b = np.zeros(self.dofs_u)
        self.load_func = load_func

        self.assembler(mesh)


    def ravthom_coeffs(self, mesh, ind_el):
        """
        Function to get the coefficients for 
        the Raviart-Thomas basis functions
        for an element. Each basis function is of the
        form
        phi[i](x) = a[i]*(x-points[element[i]]),
        i.e, a scaled point source at each vertex,
        and ravthom_coeffs(...) returns a as a list.
        """
        a = np.zeros(3)
        # iterate over edges:
        element = mesh.elements[ind_el]
        for i in range(3):
            # Get vertex index:
            ind_v = element[i]
            # Get edge index:
            ind_ed = mesh.element_edges[ind_el][i]
            # Get vertex index of the first point on the edge:
            ind_ved = mesh.faces[ind_ed][0]

            #Get the relevant normal vector:
            normal = mesh.normals[ind_ed]

            # Get inverse of coefficient:
            coeffinv = np.inner(self.points[ind_ved], normal)-\
                       np.inner(self.points[ind_v], normal)
            
            # Append to a:
            a[i] = (1./coeffinv)

        # And finally return array:
        return a
        

    def local_mass_sig(self, local_points, rt_coeffs):
        """
        Function for returning the local mass matrix for the Raviart-Thomas
        space over an element.
        """
        # Initialize local matrix
        A_loc = np.zeros( (3,3) )
        

        # Iterate over rows:
        for alpha in xrange(3):
            # Iterate over columns:
            for beta in xrange(3):
                A_loc[alpha, beta] += gauss.gaussian_quad_2d(local_points[0],\
                                                            local_points[1],\
                                                            local_points[2],\
                                                            4,\
                                           lambda x: np.inner(rt_coeffs[alpha]*(x-local_points[alpha]),\
                                                              rt_coeffs[beta]*(x-local_points[beta])))


        return A_loc

    def local_div_matrix(self, local_points, rt_coeffs):
        """
        Function for return the local divergence matrix over an element.
        """
        B_loc = np.zeros( (1, 3) )
        # Area of triangle:
        area= gauss.gaussian_quad_2d(local_points[0],\
                                     local_points[1],\
                                     local_points[2],\
                                     1,\
                                     lambda x: 1.)
        B_loc[0] = area*rt_coeffs
        return B_loc

    def local_loading(self, local_points):
        """
        Function for returning local loading function.
        Recall that we test against piecewise constants
        on one of the equations.
        """
        return -gauss.gaussian_quad_2d(local_points[0],\
                                      local_points[1],\
                                      local_points[2],\
                                      4,\
                                      self.load_func)

    ##    THE ASSEMBLER   ##
    def assembler(self, mesh):
        # We start be iterating over each element:
        for ind_el, element in enumerate(mesh.elements):
            # Get edge indices:
            edges = mesh.element_edges[ind_el]
            # Get the RT-coefficients:
            rt_coeffs = self.ravthom_coeffs(mesh,ind_el)

            # Get local points:
            local_points = np.array( self.points[element] )

            # Add contributions to mass:
            self.A[np.ix_(edges,edges)] += self.local_mass_sig(local_points, rt_coeffs)
            
            # Add contribution to divergence matrix:
            self.B[np.ix_([ind_el], edges)] += self.local_div_matrix(local_points, rt_coeffs)
            
            # Add contribution to loading:
            self.b[ind_el] += self.local_loading(local_points) 
