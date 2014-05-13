"""Module provides the functions to generate the test model for the FFF scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian


def nondimensional_parameters():
   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   return [True, True, True]


def all_fields():
   return ["streamfunction", "velocityz", "temperature"]


def implicit_fields(field_row):
   return [("streamfunction",""), ("velocityz",""), ("temperature","")]


def explicit_fields(field_row):
   return []


def equation_info(res, field_row):
   """Provide description of the system of equation"""

   # Matrix operator is real
   is_complex = False
   # Implicit field coupling
   im_fields = implicit_fields(field_row)
   # Additional explicit linear fields
   ex_fields = explicit_fields(field_row)
   # Equation doesn't have geometric coupling
   has_geometric_coupling = False
   # Rows per equation block and number of rhs
   block_info = (1, 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, block_info)


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   print("CALLING QI OPERATOR")
   print(res)
   print(eigs)
   print(bcs)
   print(field_row)

   if field_row == ("streamfunction",""):
      mat = spsp.identity(1)

   elif field_row == ("velocityz",""):
      mat = spsp.identity(1)

   elif field_row == ("temperature",""):
      mat = spsp.identity(1)

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block of linear operator"""

   if field_row == ("streamfunction",""):
      if field_row == ("streamfunction",""):
         mat = cartesian.lapl(eigs[0],eigs[1],eigs[2])

      elif field_row == ("velocityz",""):
         mat = cartesian.zblk(1)

      elif field_row == ("temperature",""):
         mat = cartesian.zblk(1)

   elif field_row == ("velocityz",""):
      if field_row == ("streamfunction",""):
         mat = cartesian.zblk(1)

      elif field_row == ("velocityz",""):
         mat = cartesian.lapl2(eigs[0],eigs[1],eigs[2])

      elif field_row == ("temperature",""):
         mat = cartesian.zblk(1)

   elif field_row == ("temperature",""):
      if field_row == ("streamfunction",""):
         mat = cartesian.zblk(1)

      elif field_row == ("velocityz",""):
         mat = cartesian.zblk(1)

      elif field_row == ("temperature",""):
         mat = cartesian.lapl(eigs[0],eigs[1],eigs[2])

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   if field_row == ("streamfunction",""):
      mat = cartesian.id()

   elif field_row == ("velocityz",""):
      mat = cartesian.lapl(eigs[0],eigs[1],eigs[2])

   elif field_row == ("temperature",""):
      mat = cartesian.id()

   return mat


def time(res, eq_params, eigs, bcs, field_row):
   """Create the time derivative operator"""

   print("CALLING TIME OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)
   
   return utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   print("CALLING IMPLICIT OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)

   return utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   print("CALLING EXPLICIT OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)
   print(field_col)

   return -linear_block(res, eq_params, eigs, field_row, field_col)
