"""Module provides the functions to generate the test model for the SLF (shell) scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets
import shell


def nondimensional_parameters():
   return ["prandtl", "rayleigh"]


def periodicity():
   return [False, False, False]


def all_fields():
   return ["velocity","temperature"]


def implicit_fields(field_row):
   return [("velocity","tor"), ("velocity","pol"), ("temperature","")]


def explicit_fields(field_row):
   return []


def equation_info(res, field_row):
   """Provide description of the system of equation"""

   is_complex = False
   im_fields = implicit_fields(field_row)
   ex_fields = explicit_fields(field_row)

   return (is_complex,im_fields,ex_fields)


def qi(res, eigs, bcs, field_row):
   """Create the quasi-inverse operator"""

   print("CALLING QI OPERATOR")
   print(res)
   print(eigs)
   print(bcs)
   print(field_row)

   if field_row == ("velocity","tor"):
      mat = shell.i2x2(res[0],eigs[0],eigs[1])

   elif field_row == ("velocity","pol"):
      mat = shell.i4x4(res[0],eigs[0],eigs[1])

   elif field_row == ("temperature",""):
      mat = shell.i2x2(res[0],eigs[0],eigs[1])

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block of linear operator"""

   if field_row == ("velocity","tor"):
      if field_col == ("velocity","tor"):
         mat = shell.i2x2lapl(res[0],eigs[0],eigs[1])

      elif field_col == ("velocity","pol"):
         mat = shell.zblk(res[0],eigs[0])

      elif field_col == ("temperature",""):
         mat = shell.zblk(res[0],eigs[0])

   elif field_row == ("velocity","pol"):
      if field_col == ("velocity","tor"):
         mat = shell.zblk(res[0],eigs[0])

      elif field_col == ("velocity","pol"):
         mat = shell.i4x4lapl2(res[0],eigs[0],eigs[1])

      elif field_col == ("temperature",""):
         mat = shell.zblk(res[0],eigs[0])

   elif field_row == ("temperature",""):
      if field_col == ("velocity","tor"):
         mat = shell.zblk(res[0],eigs[0])

      elif field_col == ("velocity","pol"):
         mat = shell.zblk(res[0],eigs[0])

      elif field_col == ("temperature",""):
         mat = shell.i2x2lapl(res[0],eigs[0],eigs[1])

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):
   """Create matrix block of time operator"""

   if field_row == ("velocity","tor"):
      mat = shell.i2x2(res[0],eigs[0],eigs[1])

   elif field_row == ("velocity","pol"):
      mat = shell.i4x4lapl(res[0],eigs[0],eigs[1])

   elif field_row == ("temperature",""):
      mat = shell.i2x2(res[0],eigs[0],eigs[1])

   return mat


def time(res, eq_params, eigs, bcs, field_row):
   """Create the time derivative operator"""

   print("CALLING TIME OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)
   
   return utils.build_diag_matrix(implicit_fields(field_row), time_block, (res,eq_params,eigs,bcs))


def implicit_linear(res, eq_params, eigs, bcs, field_row):
   """Create the implicit linear operator"""

   print("CALLING IMPLICIT OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)

   return utils.build_block_matrix(implicit_fields(field_row), linear_block, (res,eq_params,eigs,bcs))


def explicit_linear(res, eq_params, eigs, bcs, field_row, field_col):
   """Create the explicit linear operator"""

   print("CALLING EXPLICIT OPERATOR")
   print(res)
   print(eq_params)
   print(eigs)
   print(bcs)
   print(field_row)
   print(field_col)

   return -linear_block(res, eq_params, eigs, field_row, field_col)
