"""Module provides the functions to generate the Boussinesq Beta 3DQG model"""

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

   if field_row == ("vorticity",""):
      fields = []
   else:
      fields = [("streamfunction",""), ("velocityz",""), ("temperature","")]

   return fields


def explicit_fields(field_row):
   """Get the list of fields with explicit linear dependence"""

   if field_row == ("vorticity",""):
      fields = [("vorticity","")]
   else:
      fields = []

   return fields


def equation_info(res, field_row):
   """Provide description of the system of equation"""

   # Matrix operator is complex except for vorticity
   if field_row == ("vorticity",""):
      is_complex = False
   else:
      is_complex = True

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

   return (is_complex, im_fields, ex_fields, has_gometric_coupling, index_mode, block_info)


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
         chi = eq_params['chi']
         G = eq_params['gamma']
         k = eigs[0]/2

         bc = None
         bcId = bcs.get(field_col[0], -1)
         if bcId == 0:
            bc_field = {}
            bc_field[("streamfunction","")] = {'x':[40],'z':[10, -1j*k*np.tan(chi*np.pi/180)/G]}
            bc_field[("velocityz","")] = {'x':[20],'z':[11]}
            bc_field[("temperature","")] = {'x':[20],'z':[0]}
            if field_col == field_row:
               bc = bc_field[field_col]
         elif bcId == 1:
            bc_field = {}
            bc_field[("streamfunction","")] = {'x':[41],'z':[10, -1j*k*np.tan(chi*np.pi/180)/G]}
            bc_field[("velocityz","")] = {'x':[21],'z':[11]}
            bc_field[("temperature","")] = {'x':[20],'z':[0]}
            if field_col == field_row:
               bc = bc_field[field_col]

         if field_row == ("streamfunction","") and field_col == ("velocityz",""):
            bc = {'x':[0],'z':[10]}
         elif field_row == ("velocityz","") and field_col == ("streamfunction",""):
            bc = {'x':[0],'z':[11, 1j*k*np.tan(chi*np.pi/180)/G]}

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
      mat = c2d.i4j1(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("velocityz",""):
      mat = c2d.i2j1(res[0],res[2], {'x':[0], 'z':[0]})

   elif field_row == ("temperature",""):
      mat = c2d.i2j0(res[0],res[2], {'x':[0], 'z':[0]})

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col, linearize = False):
   """Create matrix block of linear operator"""

   Pr = eq_params['prandtl']
   Ra = eq_params['rayleigh']
   G = eq_params['gamma']
   k = eigs[0]/2

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      if field_col == ("streamfunction",""):
         mat = c2d.i4j1lapl2h(res[0],res[2],k, bc)

      elif field_col == ("velocityz",""):
         mat = c2d.i4j1d0d1(res[0],res[2], bc)

      elif field_col == ("temperature",""):
         mat = c2d.i4j1(res[0],res[2], bc, 1j*k*(Ra/(16*Pr)))

   elif field_row == ("velocityz",""):
      if field_col == ("streamfunction",""):
         mat = c2d.i2j1d0d1(res[0],res[2], bc, (-1/G**2))

      elif field_col == ("velocityz",""):
         mat = c2d.i2j1laplh(res[0],res[2],k, bc)

      elif field_col == ("temperature",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

   elif field_row == ("temperature",""):
      if field_col == ("streamfunction",""):
         mat = c2d.i2j0(res[0],res[2], bc, 1j*k)

      elif field_col == ("velocityz",""):
         mat = c2d.zblk(res[0],res[2],2,2, bc)

      elif field_col == ("temperature",""):
         mat = c2d.i2j0laplh(res[0],res[2],k, bc, (1/Pr))

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   k = eigs[0]/2

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_row)
   if field_row == ("streamfunction",""):
      mat = c2d.i4j1laplh(res[0],res[2],k, bc)

   elif field_row == ("velocityz",""):
      mat = c2d.i2j1(res[0],res[2], bc)

   elif field_row == ("temperature",""):
      mat = c2d.i2j0(res[0],res[2], bc)

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

   mat = -linear_block(res, eq_params, eigs, field_row, field_col)
   return mat
