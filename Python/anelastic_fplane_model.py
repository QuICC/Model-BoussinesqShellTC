"""Module provides the functions to generate the Anelastic F-Plane model"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_1d as c1d
import numpy as np


def nondimensional_parameters():
   """Get the list of nondimensional parameters"""

   return ["prandtl", "rayleigh", "taylor", "theta"]


def periodicity():
   """Get the domain periodicity"""

   return [False, True, True]


def all_fields():
   """Get the list of fields that need a configuration entry"""

   return ["velocityx","velocityy","velocityz","density","entropy","pressure","temperature"]


def implicit_fields(field_row):
   """Get the list of coupled fields in solve"""

   return [("streamfunction",""), ("velocityz",""), ("temperature","")]


def explicit_fields(field_row):
   """Get the list of fields with explicit linear dependence"""

   fields = []

   return fields


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
   index_mode = 1

   # Rows per equation block and number of rhs
   block_info = (res[0], 1)

   return (is_complex, im_fields, ex_fields, has_geometric_coupling, index_mode, block_info)


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
         bcId = bcs.get(field_col[0], -1)
         if bcId == 0:
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


def time(res, eq_params, eigs, bcs, fields):
   """Create the time derivative operator"""

   return utils.build_diag_matrix(fields, time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs, fields):
   """Create the implicit linear operator"""

   return utils.build_block_matrix(fields, linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   return -linear_block(res, eq_params, eigs, field_row, field_col)
