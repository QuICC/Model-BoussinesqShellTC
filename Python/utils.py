"""Module provides generic functions for the sparse chebyshev representation."""

import scipy.sparse as spsp

def build_diagonals(ns, nzrow, ds, offsets):
   """Build diagonals from function list and offsets"""

   diags = [0]*len(ds)
   for d,f in enumerate(ds):
      lb = max(0, -offsets[d])
      ub = min(len(ns),len(ns)-offsets[d])
      diags[d] = [f(n) if n > nzrow else 0 for n in ns[lb:ub]]

   return diags

def build_block_matrix(fields, func, func_args):

   tmp = [[0]*len(fields)]*len(fields)
   for r,field_row in enumerate(fields):
      for c,field_col in enumerate(fields):
         args = func_args + (field_row,field_col)
         tmp[r][c] = func(*args)

   return spsp.bmat(tmp)

def build_diag_matrix(fields, func, func_args):

   tmp = [0]*len(fields)
   for r,field_row in enumerate(fields):
      args = func_args + (field_row,)
      tmp[r] = func(*args)
   
   return spsp.block_diag(tmp)

def triplets(mat):
   mat = mat.tocoo();

   return list(zip(mat.row,mat.col,mat.data))
