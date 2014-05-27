"""Module provides the functions to generate the Boussinesq F-Plane 3DQG model"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_1d as c1d
import numpy as np


def nondimensional_parameters():
   """Get the list of nondimensional parameters"""

   return ["prandtl", "rayleigh", "theta"]


def periodicity():
   """Get the domain periodicity"""

   return [False, True, True]


def all_fields():
   """Get the list of fields that need a configuration entry"""

   return ["streamfunction", "velocityz", "temperature"]


def implicit_fields(field_row):
   """Get the list of coupled fields in solve"""

   if field_row == ("vorticity",""):
      fields = [("vorticity","")]

   elif field_row == ("meantemperature",""):
      fields = [("meantemperature","")]

   else:
      fields =  [("streamfunction",""), ("velocityz",""), ("temperature","")]

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

   # Matrix operator is complex except for vorticity and mean temperature
   if field_row == ("vorticity","") or field_row == ("meantemperature",""):
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
   index_mode = 1

   # Rows per equation block and number of rhs
   block_info = (res[0], 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, index_mode, block_info)


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
         eta2 = np.sin(np.pi*eq_params['theta']/180)
         eta3 = np.cos(np.pi*eq_params['theta']/180)
         kx = eigs[0]
         ky = eigs[1]

         bc = None
         bcId = bcs.get(field_col[0], -1)
         if bcId == 0:
            bc_field = {}
            bc_field[("streamfunction","")] = [10, -1j*kx*eta2]
            bc_field[("velocityz","")] = [11, eta3]
            bc_field[("temperature","")] = [0]
            if field_col == field_row:
               bc = bc_field[field_col]

         if field_row == ("streamfunction","") and field_col == ("velocityz",""):
            bc = [10, eta3]
         elif field_row == ("velocityz","") and field_col == ("streamfunction",""):
            bc = [11, -1j*kx*eta2]
         elif field_row == ("meantemperature","") and field_col == field_row:
            bc = [10]

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
      mat = c1d.i4(res[0], [0])

   elif field_row == ("velocityz",""):
      mat = c1d.i4(res[0], [0])

   elif field_row == ("temperature",""):
      mat = c1d.i2(res[0], [0])

   elif field_row == ("meantemperature",""):
      if eigs[0] == 0 and eigs[1] == 0:
         stat = c1d.zblk(res[0],0, [0])
         stat[0,::2] = [n/(n**2-1) - 1/(n-1) for n in np.arange(0,res[0],2)]
         tmp = c1d.zblk(res[0],0, [0])
         tmp[::2,0] = 1
         mat = c1d.i2(res[0], [0])*(c1d.qid(res[0],0,[0]) - tmp*stat)
      else:
         mat = c1d.zblk(res[0],0, [0])

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block linear operator"""

   Pr = eq_params['prandtl']
   Ra = eq_params['rayleigh']
   eta2 = np.sin(np.pi*eq_params['theta']/180)
   eta3 = np.cos(np.pi*eq_params['theta']/180)
   kx = eigs[0]
   ky = eigs[1]

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_col)
   if field_row == ("streamfunction",""):
      if field_col == ("streamfunction",""):
         mat = c1d.i1(res[0], bc, (kx**2 + (1/eta3**2)*ky**2)**2)

      elif field_col == ("velocityz",""):
         mat = c1d.i1d1(res[0], bc, 2*eta3)

      elif field_col == ("temperature",""):
         mat = c1d.i1(res[0], bc, -1j*kx*eta2*(Ra/Pr))

   elif field_row == ("velocityz",""):
      if field_col == ("streamfunction",""):
         mat = c1d.i1d1(res[0], bc, -2*eta3)

      elif field_col == ("velocityz",""):
         mat = c1d.i1(res[0], bc, -(kx**2 + (1/eta3**2)*ky**2))

      elif field_col == ("temperature",""):
         mat = c1d.i1(res[0], bc, eta3*Ra/Pr)

   elif field_row == ("temperature",""):
      if field_col == ("streamfunction",""):
         #mat = c1d.zblk(res[0],0, bc)
         mat = c1d.qid(res[0],0, bc, -1j*kx*eta2)

      elif field_col == ("velocityz",""):
         #mat = c1d.zblk(res[0],0, bc)
         mat = c1d.qid(res[0],0, bc, eta3)

      elif field_col == ("temperature",""):
         mat = c1d.qid(res[0],0, bc, -(1/Pr)*(kx**2 + (1/eta3**2)*ky**2))

   elif field_row == ("vorticity",""): 
      if field_col == field_row:
         mat = c1d.qid(res[0],0, bc, -(kx**2 + (1/eta3**2)*ky**2))

      else:
         mat = c1d.zblk(res[0],0, bc)

   elif field_row == ("meantemperature",""): 
      if field_col == field_row:
         mat = c1d.i1d1(res[0], bc, 2.0)

      else:
         mat = c1d.zblk(res[0],1, bc)

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   eta2 = np.sin(np.pi*eq_params['theta']/180)
   eta3 = np.cos(np.pi*eq_params['theta']/180)
   kx = eigs[0]
   ky = eigs[1]

   bc = convert_bc(eq_params,eigs,bcs,field_row,field_row)
   if field_row == ("streamfunction",""):
      mat = c1d.i1(res[0], bc, -(kx**2 + (1/eta3**2)*ky**2))

   elif field_row == ("velocityz",""):
      mat = c1d.i1(res[0], bc)

   elif field_row == ("temperature",""):
      mat = c1d.qid(res[0],0, bc)

   return mat


def time(res, eq_params, eigs, bcs, field_row):
   """Create the time derivative operator"""

   mat = utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))
   return mat


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   mat = utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))
   return mat


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   mat = -linear_block(res, eq_params, eigs, field_row, field_col)
   return mat
