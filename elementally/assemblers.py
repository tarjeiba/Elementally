import numpy as np
import numpy.linalg as la

from integrators import gaussian as gauss
from integrators import gl_quad_functions as gl
import meshers

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


        # Inverting the matrix of coordinates results in the coefficients for
        # the three test functions that are non-zero on this element

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
            local_points = np.array(self.points[element])
            coords = np.ones( (3,3) )
            coords[:,1:] = local_points
            coeffs = la.inv(coords)
            
            self.A[np.ix_(element, element)] += self.local_stiffness(local_points, coeffs)
            self.M[np.ix_(element, element)] += self.local_mass(local_points, coeffs)
            self.b[element] += self.local_loading(local_points, self.load_func, coeffs)


