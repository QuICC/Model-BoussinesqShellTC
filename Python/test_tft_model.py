"""Module provides the functions to generate the Boussinesq rotating convection model matrices in a sphere"""

import scipy.sparse as spsp
import utils
from utils import triplets

all_fields = [("streamfunction","scalar"), ("velocityz","scalar"), ("temperature","scalar")]
implicit_fields = [("streamfunction","scalar"), ("velocityz","scalar"), ("temperature","scalar")]
explicit_fields = []

def qi(res, eigs, bcs, field):
   """Create the quasi-inverse operator"""

   if field[0] == 'velocity' and field[1] == 'toroidal':
      mat = spsp.identity(res[0])

   elif field[0] == 'velocity' and field[1] == 'poloidal':
      mat = spsp.identity(res[0])

   elif field[0] == 'temperature' and field[1] == 'scalar':
      mat = spsp.identity(res[0])

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block for field"""

   mat = spsp.identity(res[0])

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):

   mat = spsp.identity(res[0])

   return mat


def time(res, eq_params, eigs, bcs):
   """Create the time derivative operator"""
   
   return utils.build_diag_matrix(implicit_fields, time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs):
   """Create the implicit linear operator"""

   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)

   return utils.build_block_matrix(implicit_fields, linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   return -linear_block(res, eq_params, eigs, field_row, field_col)
