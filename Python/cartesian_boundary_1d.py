"""Module provides functions to generate the boundary conditions in a cartesian 1D geometry"""

from __future__ import division

import numpy as np

def constrain(mat, bc, eq_zrows):
   """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

   if bc[0] > 0:
      bc_mat = apply_tau(mat, bc)
   elif bc[0] < 0:
      bc_mat = apply_galerkin(mat, bc, eq_zrows)
   else:
      bc_mat = mat

   return bc_mat

def apply_tau(mat, bc):
   """Add Tau lines to the matrix"""

   if bc[0] == 10:
      cond = tau_no_slip(mat.shape[0], 1)
   elif bc[0] == 11:
      cond = tau_no_slip(mat.shape[0], -1)
   elif bc[0] == 12:
      cond = tau_stress_free(mat.shape[0], 1)
   elif bc[0] == 13:
      cond = tau_stress_free(mat.shape[0], -1)
   elif bc[0] == 20:
      cond = tau_no_slip(mat.shape[0], 0)
   elif bc[0] == 21:
      cond = tau_stress_free(mat.shape[0], 0)
   elif bc[0] == 40:
      cond = tau_no_pen_ns(mat.shape[0], 0)

   bc_mat = mat
   bc_mat[0:cond.shape[0],:] = cond

   return bc_mat

def tau_no_slip(nx, pos):
   """Create the no-slip tau line(s)"""

   cond = []
   if pos >= 0:
      cond.append([1.0]*nx)

   if pos <= 0:
      cond.append([(-1.0)**i for i in np.arange(0,nx)])

   return np.array(cond)

def tau_stress_free(nx, pos):
   """Create the stress-free tau line(s)"""

   cond = []
   if pos >= 0:
      cond.append([i**2 for i in np.arange(0,nx)])

   if pos <= 0:
      cond.append([(-1.0)**(i+1)*i**2 for i in np.arange(0,nx)])

   return np.array(cond)

def tau_no_pen_ns(nx, pos):
   """Create the no penetration and no-slip tau line(s)"""

   cond = []
   if pos >= 0:
      cond.append([1.0]*nx)
      cond.append([(1/3)*(i**4 - i**2) for i in np.arange(0,nx)])

   if pos <= 0:
      cond.append([(-1.0)**i for i in np.arange(0,nx)])
      cond.append([((-1.0)**i/3)*(i**4 - i**2) for i in np.arange(0,nx)])

   return np.array(cond)

def apply_galerkin(mat, bc, eq_zero_rows):
   """Apply a Galerkin stencil on the matrix"""

   return mat
