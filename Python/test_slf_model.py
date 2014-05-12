"""Module provides the functions to generate the test model for the SLF (shell) scheme"""

import scipy.sparse as spsp
import utils
from utils import triplets

def all_fields():
   return [("streamfunction",""), ("velocityz",""), ("temperature","")]


def implicit_fields(field_row):
   return [("streamfunction",""), ("velocityz",""), ("temperature","")]


def explicit_fields(field_row):
   return []

def equation_info(res, field_row):
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

   if field_row == ("streamfunction",""):
      mat = spsp.identity(res[0]*res[2])

   elif field_row == ("velocityz",""):
      mat = spsp.identity(res[0]*res[2])

   elif field_row == ("temperature",""):
      mat = spsp.identity(res[0]*res[2])

   return mat


def linear_block(res, eq_params, eigs, bcs, field_row, field_col):
   """Create matrix block for field"""

   mat = spsp.identity(res[0]*res[2])

   return mat


def time_block(res, eq_params, eigs, bcs, field_row):

   mat = spsp.identity(res[0]*res[2])

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
