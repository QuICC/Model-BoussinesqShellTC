"""Module provides the functions to generate the test model for the TFF scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_1d as c1d


def nondimensional_parameters():
   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   return [False, True, True]


def all_fields():
   return ["streamfunction","velocityz""temperature"]


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
   block_info = (res[0], 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, block_info)


def convert_bc(eq_params, eigs, bcs, field_row, field_col):
   """Convert simulation input boundary conditions to ID"""

   if field_col == ("streamfunction",""):
      if bcs[field_col] == 0:
         if field_row == ("streamfunction",""):
            bc = [40]
         elif field_row == ("velocityz",""):
            bc = [0]
         elif field_row == ("temperature",""):
            bc = [0]
   elif field_col == ("velocityz",""):
      if bcs[field_col] == 0:
         if field_row == ("streamfunction",""):
            bc = [0]
         elif field_row == ("velocityz",""):
            bc = [20]
         elif field_row == ("temperature",""):
            bc = [0]
   elif field_col == ("temperature",""):
      if bcs[field_col] == 0:
         if field_row == ("streamfunction",""):
            bc = [0]
         elif field_row == ("velocityz",""):
            bc = [0]
         elif field_row == ("temperature",""):
            bc = [20]
   
   return bc


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   print("CALLING QI OPERATOR")
   print(res)
   print(eigs)
   print(bcs)
   print(field_row)

   if field_row == ("streamfunction",""):
      mat = c1d.i2(res[0], [0])

   elif field_row == ("velocityz",""):
      mat = c1d.i4(res[0], [0])

   elif field_row == ("temperature",""):
      mat = c1d.i2(res[0], [0])

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block linear operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      if field_row == ("streamfunction",""):
         mat = c1d.i4lapl2(res[0],eigs[0],eigs[1], bc)

      elif field_row == ("velocityz",""):
         mat = c1d.zblk(res[0], bc)

      elif field_row == ("temperature",""):
         mat = c1d.zblk(res[0], bc)

   elif field_row == ("velocityz",""):
      if field_row == ("streamfunction",""):
         mat = c1d.zblk(res[0], bc)

      elif field_row == ("velocityz",""):
         mat = c1d.i2lapl(res[0],eigs[0],eigs[1], bc)

      elif field_row == ("temperature",""):
         mat = c1d.zblk(res[0] ,bc)

   elif field_row == ("temperature",""):
      if field_row == ("streamfunction",""):
         mat = c1d.zblk(res[0], bc)

      elif field_row == ("velocityz",""):
         mat = c1d.zblk(res[0], bc)

      elif field_row == ("temperature",""):
         mat = c1d.i2lapl(res[0],eigs[0],eigs[1], bc)

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_row)
   if field_row == ("streamfunction",""):
      mat = c1d.i4lapl(res[0],eigs[0],eigs[1], bc)

   elif field_row == ("velocityz",""):
      mat = c1d.i2(res[0], bc)

   elif field_row == ("temperature",""):
      mat = c1d.i2(res[0], bc)

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
