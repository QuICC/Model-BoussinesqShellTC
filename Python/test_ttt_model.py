"""Module provides the functions to generate the test model for the TTT scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_3d as c3d


def nondimensional_parameters():
   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   return [False, False, False]


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
   # Index mode: SLOWEST = 0, MODE = 1
   index_mode = 0

   # Rows per equation block and number of rhs
   block_info = (res[0]*res[1]*res[2], 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, block_info)


def convert_bc(eq_params, eigs, bcs, field_row, field_col):
   """Convert simulation input boundary conditions to ID"""

   use_tau_boundary = True
   # Impose no boundary conditions
   no_bc = {'x':[0],'y':[0],'z':[0]}
   if bcs["bcType"] == 2:
      bc = no_bc
   else:
      # Impose no boundary conditions
      if bcs["bcType"] == 1 and use_tau_boundary:
         bc = no_bc
      else: #bcType == 0 or Galerkin boundary
         bc = None
         if bcs[field_col[0]] == 0:
            bc_field = {}
            bc_field[("streamfunction","")] = {'x':[40],'y':[40],'z':[40]}
            bc_field[("velocityz","")] = {'x':[20],'y':[20],'z':[20]}
            bc_field[("temperature","")] = {'x':[20],'y':[20],'z':[20]}
            if field_col == field_row:
               bc = bc_field[field_col]

         if bc is None:
            if use_tau_boundary:
               bc = no_bc
            else:
               bc = {}
               for k,v in bc_field[field_col]:
                  bc[k] = v
                  bc[k][0] = -v[0]
   
   return bc


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   print("CALLING QI OPERATOR")

   if field_row == ("streamfunction",""):
      mat = c3d.i2j2k2(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

   elif field_row == ("velocityz",""):
      mat = c3d.i4j4k4(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

   elif field_row == ("temperature",""):
      mat = c3d.i2j2k2(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block of linear operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      if field_row == ("streamfunction",""):
         mat = c3d.i2j2k2lapl(res[0],res[1],res[2], bc)

      elif field_row == ("velocityz",""):
         mat = c3d.zblk(res[0],res[1],res[2],4,4,4, bc)

      elif field_row == ("temperature",""):
         mat = c3d.zblk(res[0],res[1],res[2],4,4,4, bc)

   elif field_row == ("velocityz",""):
      if field_row == ("streamfunction",""):
         mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

      elif field_row == ("velocityz",""):
         mat = c3d.i4j4k4lapl2(res[0],res[1],res[2], bc)

      elif field_row == ("temperature",""):
         mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

   elif field_row == ("temperature",""):
      if field_row == ("streamfunction",""):
         mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

      elif field_row == ("velocityz",""):
         mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

      elif field_row == ("temperature",""):
         mat = c3d.i2j2k2lapl(res[0],res[1],res[2], bc)

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      mat = c3d.i2j2k2(res[0],res[1],res[2], bc)

   elif field_row == ("velocityz",""):
      mat = c3d.i4j4k4lapl(res[0],res[1],res[2], bc)

   elif field_row == ("temperature",""):
      mat = c3d.i2j2k2(res[0],res[1],res[2], bc)

   return mat


def time(res, eq_params, eigs, bcs, field_row):
   """Create the time derivative operator"""

   print("CALLING TIME OPERATOR")
   
   return utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   print("CALLING IMPLICIT OPERATOR")

   return utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   print("CALLING EXPLICIT OPERATOR")

   return -linear_block(res, eq_params, eigs, field_row, field_col)
