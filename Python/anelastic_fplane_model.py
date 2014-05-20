"""Module provides the functions to generate the Anelastic F-Plane model"""

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

   use_tau_boundary = True
   # Impose no boundary conditions
   no_bc = [0]
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
            bc_field[("streamfunction","")] = [40]
            bc_field[("velocityz","")] = [20]
            bc_field[("temperature","")] = [20]
            if field_col == field_row:
               bc = bc_field[field_col]

         if bc is None:
            if use_tau_boundary:
               bc = no_bc
            else:
               bc = -bc_field[field_col]
   
   return bc


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

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
         mat = c1d.zblk(res[0],4, bc)

      elif field_row == ("temperature",""):
         mat = c1d.zblk(res[0],4, bc)

   elif field_row == ("velocityz",""):
      if field_row == ("streamfunction",""):
         mat = c1d.zblk(res[0],2, bc)

      elif field_row == ("velocityz",""):
         mat = c1d.i2lapl(res[0],eigs[0],eigs[1], bc)

      elif field_row == ("temperature",""):
         mat = c1d.zblk(res[0],2, bc)

   elif field_row == ("temperature",""):
      if field_row == ("streamfunction",""):
         mat = c1d.zblk(res[0],2, bc)

      elif field_row == ("velocityz",""):
         mat = c1d.zblk(res[0],2, bc)

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

   return utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   return utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   return -linear_block(res, eq_params, eigs, field_row, field_col)
