"""Module provides the functions to generate the test model for the TFT scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_2d as c2d


def nondimensional_parameters():
   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   return [False, True, False]


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
   block_info = (res[0]*res[2], 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, block_info)


def convert_bc(eq_params, eigs, bcs, field_row, field_col):
   """Convert simulation input boundary conditions to ID"""

   if field_col == ("streamfunction",""):
      if bcs[field_col[0]] == 0:
         if field_row == ("streamfunction",""):
            bc = {'x':[40],'z':[40]}
         elif field_row == ("velocityz",""):
            bc = {'x':[0],'z':[0]}
         elif field_row == ("temperature",""):
            bc = {'x':[0],'z':[0]}
   elif field_col == ("velocityz",""):
      if bcs[field_col[0]] == 0:
         if field_row == ("streamfunction",""):
            bc = {'x':[0],'z':[0]}
         elif field_row == ("velocityz",""):
            bc = {'x':[20],'z':[20]}
         elif field_row == ("temperature",""):
            bc = {'x':[0],'z':[0]}
   elif field_col == ("temperature",""):
      if bcs[field_col[0]] == 0:
         if field_row == ("streamfunction",""):
            bc = {'x':[0],'z':[0]}
         elif field_row == ("velocityz",""):
            bc = {'x':[0],'z':[0]}
         elif field_row == ("temperature",""):
            bc = {'x':[20],'z':[20]}
   
   return bc


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   print("CALLING QI OPERATOR")

   if field_row == ("streamfunction",""):
      mat = c2d.i4j4(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("velocityz",""):
      mat = c2d.i2j2(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("temperature",""):
      mat = c2d.i2j2(res[0],res[2], {'x':[0], 'z':[0]})

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block of linear operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      if field_col == ("streamfunction",""):
         mat = c2d.i4j4lapl2(res[0],res[2],eigs[0], bc)

      elif field_col == ("velocityz",""):
         mat = c2d.zblk(res[0],res[2],4,4, bc)

      elif field_col == ("temperature",""):
         mat = c2d.zblk(res[0],res[2],4,4, bc)

   elif field_row == ("velocityz",""):
      if field_col == ("streamfunction",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

      elif field_col == ("velocityz",""):
         mat = c2d.i2j2lapl(res[0],res[2],eigs[0], bc)

      elif field_col == ("temperature",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

   elif field_row == ("temperature",""):
      if field_col == ("streamfunction",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

      elif field_col == ("velocityz",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

      elif field_col == ("temperature",""):
         mat = c2d.i2j2lapl(res[0],res[2],eigs[0], bc)

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_row)
   if field_row == ("streamfunction",""):
      mat = c2d.i4j4lapl(res[0],res[2],eigs[0], bc)

   elif field_row == ("velocityz",""):
      mat = c2d.i2j2(res[0],res[2], bc)

   elif field_row == ("temperature",""):
      mat = c2d.i2j2(res[0],res[2], bc)

   return mat


def time(res, eq_params, eigs, bcs, field_row):
   """Create the time derivative operator"""

   print("CALLING TIME OPERATOR")
   
   mat = utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))

   return mat


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   print("CALLING IMPLICIT OPERATOR")

   mat = utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))

   return mat


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   print("CALLING EXPLICIT OPERATOR")

   mat = -linear_block(res, eq_params, eigs, field_row, field_col)

   return mat
