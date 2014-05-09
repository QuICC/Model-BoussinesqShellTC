"""Module provides the functions to generate the Boussinesq rotating convection model matrices in a sphere"""

import numpy as np
import scipy.sparse as spsp

all_fields = [("velocity","toroidal"), ("velocity","toroidal"), ("temperature","scalar")]
implicit_fields = [("velocity","toroidal"), ("velocity","toroidal"), ("temperature","scalar")]
explicit_fields = []

def qi(res, eq_params, eigs, field, bcs):
   """Create the quasi-inverse operator"""

   if field[0] == 'velocity' and field[1] == 'toroidal':
      mat = spsp.identity(res[0])

   elif field[0] == 'velocity' and field[1] == 'poloidal':
      mat = spsp.identity(res[0])

   elif field[0] == 'temperature' and field[1] == 'scalar':
      mat = spsp.identity(res[0])

   return mat


def time(res, eq_params, eigs, bcs):
   """Create the time derivative operator"""
   
   tor_time = spsp.identity(res[0])
   pol_time = spsp.identity(res[0])
   temp_time = spsp.identity(res[0])

   return spsp.block_diag((tor_time,pol_time,temp_time))


def implicit_linear(res, eq_params, eigs, bcs):
   """Create the implicit linear operator"""

   tmp = [[0]*len(implicit_fields)]*len(implicit_fields)
   for r,field_row in enumerate(implicit_fields):
      for c,field_col in enumerate(implicit_fields):
         tmp[r][c] = block(res, eq_params, eigs, field_row, field_col)

   return spsp.bmat(tmp)


def explicit_linear(res, eq_params, eigs, field_row, field_col, bcs):
   """Create the explicit linear operator"""

   return -block(res, eq_params, eigs, field_row, field_col)


def block(res, eq_params, eigs, field_row, field_col, bcs):
   """Create matrix block for field"""

   if field_row[0] == 'velocity' and field_row[1] == 'toroidal':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = spsp.identity(res[0])

   elif field_row[0] == 'velocity' and field_row[1] == 'poloidal':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = spsp.identity(res[0])

   elif field_row[0] == 'temperature' and field_row[1] == 'scalar':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = spsp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = spsp.identity(res[0])

   return mat


