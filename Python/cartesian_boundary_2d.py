"""Module provides functions to generate the boundary conditions in a cartesian 2D geometry"""

from __future__ import division

import numpy as np
import scipy.sparse as spsp
import cartesian_boundary_1d as c1dbc

def qid(nx, q, bc):
   """Create a quasi indentity"""

   if bc[0] < 0:
      mat = spsp.identity(nx-bc[0]//10)
   else:
      offsets = [0]
      diags = [[0]*q + [1]*(nx-q)]

      mat = spsp.diags(diags, offsets)

   return mat.tocsr()


def constrain(mat, nx, nz, bc, eq_zrows_x, eq_zrows_z):
   """Contrain the matrix with the Tau boundary condition"""

   bc_mat = mat
   if bc['x'][0] > 0:
      bid = qid(nx,0,[0])
      bid = c1dbc.constrain(bid, bc['x'], 0)
      bc_mat = bc_mat + spsp.kron(qid(nz,0,bc['z']), bid)

   if bc['z'][0] > 0:
      bid = qid(nz,0,[0])
      bid = c1dbc.constrain(bid, bc['z'], 0)
      bc_mat = bc_mat + spsp.kron(bid, qid(nx,eq_zrows_x,bc['x']))
  
   return bc_mat

