"""Module provides the functions to generate the ???? model"""

import scipy.sparse as spsp
import utils
from utils import triplets
import cartesian_1d as c1d
import numpy as np


def nondimensional_parameters():
   """Get the list of nondimensional parameters"""

   return []


def periodicity():
   """Get the domain periodicity"""

   return [?, ?, ?]


def all_fields():
   """Get the list of fields that need a configuration entry"""

   return []


def implicit_fields(field_row):
   """Get the list of coupled fields in solve"""

   fields = []

   return fields


def explicit_fields(field_row):
   """Get the list of fields with explicit linear dependence"""

   fields = []

   return fields


def equation_info(res, field_row):
   """Provide description of the system of equation"""

   # Matrix operator is complex?
   is_complex = ?

   # Implicit field coupling
   im_fields = implicit_fields(field_row)
   # Additional explicit linear fields
   ex_fields = explicit_fields(field_row)

   # Equation doesn't have geometric coupling
   has_geometric_coupling = False
   # Index mode: SLOWEST = 0, MODE = 1
   index_mode = ?

   # Rows per equation block and number of rhs
   block_info = (?, ?)

   return (is_complex,im_fields,ex_fields,has_geometric_coupling, index_mode, block_info)


def convert_bc(eq_params, eigs, bcs, field_row, field_col):
   """Convert simulation input boundary conditions to ID"""

   use_tau_boundary = True

   bc = []
   
   return bc


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block linear operator"""

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

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
