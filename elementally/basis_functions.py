import numpy as np

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
