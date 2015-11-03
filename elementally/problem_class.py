import numpy as np
import numpy.linalg as la

import meshers
#################################################
##
##            BASE PROBLEM CLASS:
##
#################################################
class FunctionSpaceBase2D:
    """Base class for function spaces in 2D.
    Should contain the following:
      Mesh: Mesh of domain where PDE is to be solved.
      DOF map: A local to global mapping of degrees of freedom.
      DOFs of dim: An array containing number of DOFs per k-dimensional
        subsimplex. Let's call it dof_dim.
        For instance dof_dim[0] is number of dofs per vertex, 
        dof_dim[1] is number of dofs per face,
        dof_dim[2] is number of dofs per triangle.
      dofs_per_el: Integer, number of dofs per elements.
      funcs: Functions for each dof, on a reference triangle.
      funcders: Derivative of the above functions.
    """
    
    def __init__(self, mesh):
        self.mesh = mesh  # Supposed to be instance of ElementallyMeshInfo.

################################################
##
##          CONTINUOUS LINEAR ELEMENTS:
##
################################################
class ContinuousLinear(FunctionSpaceBase2D):
    # Constructor is inherited from ProblemBase2D.
     
    # dofmap is same as element connectivity in this case: 
    dofmap = np.array(self.mesh.elements)
  
    # dof_dim array: One per vertex, rest is zero.
    dof_dim = np.array( (1, 0, 0) )
    dofs_per_el = 3
    
    funcs = [cg_phi0, cg_phi1, cg_phi2]
    funcders = [cg_grad_phi0, \
                cg_grad_phi1, \
                cg_grad_phi2]

#################################################
##
##        RAVIART-THOMAS ELEMENTS:
##
#################################################
class RaviartThomas(FunctionSpaceBase2D):
    # Constructor inherited from ProblemBase2D.

    # dof_dim array: None per vertex, one per edge, none per triangle.
    dof_dim = np.array( (0, 1, 0) )
    dof_per_el = 3




##################################################
##
##    BASIS FUNCTIONS ON REFERENCE TRIANGLE:
##
##################################################
# For continuous linear elements:
def cg_phi0(x):
    return 1.-x[0] - x[1]

def cg_phi1(x):
    return x[0]

def cg_phi2(x):
    return x[1]

def cg_grad_phi0(x):
    return np.array((-1., -1.))

def cg_grad_phi1(x):
    return np.array((1., 0.))

def cg_grad_phi2(x):
    return np.array((0., 1.))

# For Raviart-Thomas:
def rt_phi0(x):
  return np.sqrt(2)*np.array( (x[0], y[0]) )

def rt_phi1(x):
  return np.array( (x[0]-1., x[1]) )

def rt_phi2(x):
  return np.array( (x[0], x[1]-1.) )

def rt_div_phi0(x):
  return 2.*np.sqrt(2)

def rt_div_phi1(x):
  return 2.

def rt_div_phi2(x):
  return 2.
