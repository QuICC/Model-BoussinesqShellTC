"""Module provides the functions to generate the Boussinesq rotating convection model matrices in a sphere"""

import numpy as np
import scipy.sparse as sp

all_fields = [("velocity","toroidal"), ("velocity","toroidal"), ("temperature","scalar")]
implicit_fields = [("velocity","toroidal"), ("velocity","toroidal"), ("temperature","scalar")]
explicit_fields = []

def qi(res, eq_params, eigs, field):
   """Create the quasi-inverse operator"""

   if field[0] == 'velocity' and field[1] == 'toroidal':
      mat = sp.identity(res[0])

   elif field[0] == 'velocity' and field[1] == 'poloidal':
      mat = sp.identity(res[0])

   elif field[0] == 'temperature' and field[1] == 'scalar':
      mat = sp.identity(res[0])

   return mat


def time(res, eq_params, eigs):
   """Create the time derivative operator"""
   
   tor_time = sp.identity(res[0])
   pol_time = sp.identity(res[0])
   temp_time = sp.identity(res[0])

   return sp.block_diag((tor_time,pol_time,temp_time))


def linear(res, eq_params, eigs):
   """Create the implicit linear operator"""

   tmp = [[0]*len(implicit_fields)]*len(implicit_fields)
   for r,field_row in enumerate(implicit_fields):
      for c,field_col in enumerate(implicit_fields):
         tmp[r][c] = block(res, eq_params, eigs, field_row, field_col)

   return sp.bmat(tmp)


def explicit(res, eq_params, eigs, field_row, field_col):
   """Create the explicit linear operator"""

   return -block(res, eq_params, eigs, field_row, field_col)


def block(res, eq_params, eigs, field_row, field_col):
   """Create matrix block for field"""

   if field_row[0] == 'velocity' and field_row[1] == 'toroidal':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = sp.identity(res[0])

   elif field_row[0] == 'velocity' and field_row[1] == 'poloidal':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = sp.identity(res[0])

   elif field_row[0] == 'temperature' and field_row[1] == 'scalar':
      if field_col[0] == 'velocity' and field_col[1] == 'toroidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'velocity' and field_col[1] == 'poloidal':
         mat = sp.identity(res[0])
      elif field_col[0] == 'temperature' and field_col[1] == 'scalar':
         mat = sp.identity(res[0])

   return mat


