"""Module provides the functions to generate the test model for the TFT scheme"""

from __future__ import division

import numpy as np
import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_2d as c2d
import scipy.io as io


def nondimensional_parameters():
   """Get the list of nondimensional parameters"""

   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   """Get the domain periodicity"""

   return [False, True, False]


def all_fields():
   """Get the list of fields that need a configuration entry"""

   return ["streamfunction", "velocityz", "temperature"]


def implicit_fields(field_row):
   """Get the list of coupled fields in solve"""

   # Solve as splitted equations
   if True:
      fields = [("streamfunction",""), ("velocityz",""), ("temperature","")]

   # Solve as coupled equations
   else:
      fields = [field_row]

   return fields


def explicit_fields(field_row):
   """Get the list of fields with explicit linear dependence"""

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
   block_info = (res[0]*res[2], 1)

   return (is_complex, im_fields, ex_fields, has_geometric_coupling, index_mode, block_info)


def convert_bc(eq_params, eigs, bcs, field_row, field_col):
   """Convert simulation input boundary conditions to ID"""

   use_tau_boundary = True

   # Impose no boundary conditions
   no_bc = {'x':[0],'z':[0]}
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
            bc_field[("streamfunction","")] = {'x':[40],'z':[40]}
            bc_field[("velocityz","")] = {'x':[20],'z':[20]}
            bc_field[("temperature","")] = {'x':[20],'z':[20]}
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

   if field_row == ("streamfunction",""):
      mat = c2d.i4j4(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("velocityz",""):
      mat = c2d.i2j2(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("temperature",""):
      mat = c2d.i2j2(res[0],res[2], {'x':[0], 'z':[0]})

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col, linearize = False):
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


def time(res, eq_params, eigs, bcs, fields):
   """Create the time derivative operator"""

   mat = utils.build_diag_matrix(fields, time_block, (res,eq_params,eigs,bcs))
   io.mmwrite("matrix_time_" + str(bcs["bcType"]) + "_"+ str(eigs[0]) + ".mtx", mat)
   return mat


def implicit_linear(res, eq_params, eigs, bcs, fields):
   """Create the implicit linear operator"""

   mat = utils.build_block_matrix(fields, linear_block, (res,eq_params,eigs,bcs))
   io.mmwrite("matrix_linear_" + str(bcs["bcType"]) + "_"+ str(eigs[0]) + ".mtx",mat)
   return mat


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   mat = linear_block(res, eq_params, eigs, field_row, field_col)
   return mat
