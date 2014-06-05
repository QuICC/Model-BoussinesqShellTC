"""Module provides the functions to generate the test model for the FFF scheme"""

from __future__ import division

import numpy as np
import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_0d as c0d


def nondimensional_parameters():
   """Get the list of nondimensional parameters"""

   return ["prandtl", "rayleigh", "gamma", "chi"]


def periodicity():
   """Get the domain periodicity"""

   return [True, True, True]


def all_fields():
   """Get the list of fields that need a configuration entry"""

   return ["streamfunction", "velocityz", "temperature"]


def implicit_fields(field_row):
   """Get the list of coupled fields in solve"""

   return [("streamfunction",""), ("velocityz",""), ("temperature","")]


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
   index_mode = 1

   # Rows per equation block and number of rhs
   block_info = (1, 1)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, index_mode, block_info)


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   if field_row == ("streamfunction",""):
      mat = c0d.qid()

   elif field_row == ("velocityz",""):
      mat = c0d.qid()

   elif field_row == ("temperature",""):
      mat = c0d.qid()

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col, linearize = False):
   """Create matrix block of linear operator"""

   if field_row == ("streamfunction",""):
      if field_col == ("streamfunction",""):
         mat = c0d.lapl(eigs[0],eigs[1],eigs[2])

      elif field_col == ("velocityz",""):
         mat = c0d.zblk()

      elif field_col == ("temperature",""):
         mat = c0d.zblk()

   elif field_row == ("velocityz",""):
      if field_col == ("streamfunction",""):
         mat = c0d.zblk(0)

      elif field_col == ("velocityz",""):
         mat = c0d.lapl2(eigs[0],eigs[1],eigs[2])

      elif field_col == ("temperature",""):
         mat = c0d.zblk()

   elif field_row == ("temperature",""):
      if field_col == ("streamfunction",""):
         mat = c0d.zblk()

      elif field_col == ("velocityz",""):
         mat = c0d.zblk()

      elif field_col == ("temperature",""):
         mat = c0d.lapl(eigs[0],eigs[1],eigs[2])

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   if field_row == ("streamfunction",""):
      mat = c0d.qid()

   elif field_row == ("velocityz",""):
      mat = c0d.lapl(eigs[0],eigs[1],eigs[2])

   elif field_row == ("temperature",""):
      mat = c0d.qid()

   return mat


def time(res, eq_params, eigs, bcs, fields):
   """Create the time derivative operator"""

   mat = utils.build_diag_matrix(fields, time_block, (res,eq_params,eigs,bcs))
   return mat


def implicit_linear(res, eq_params, eigs, bcs, fields):
   """Create the implicit linear operator"""

   mat = utils.build_block_matrix(fields, linear_block, (res,eq_params,eigs,bcs))
   return mat


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   mat = -linear_block(res, eq_params, eigs, field_row, field_col)
   return mat
