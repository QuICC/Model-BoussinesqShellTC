"""Module provides functions to generate sparse operators in a cartesian box with two periodic directions."""

from __future__ import division

import numpy as np
import scipy.sparse as spsp
import utils
import cartesian_boundary_1d as c1dbc


def zblk(nx, q, bc):
   """Create a block of zeros"""

   mat = spsp.coo_matrix((nx,nx))
   return c1dbc.constrain(mat,bc,q)


def i2d2(nx, bc, coeff = 1.0):
   """Create a quasi identity block of order 2"""

   return qid(nx,2, bc, coeff)


def i1(nx, bc, coeff = 1.0):
   """Create operator for 1st integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-1,2,2)
   nzrow = 0

   # Generate 1st subdiagonal
   def d_1(n):
      return 1.0/(2*n)

   # Generate 1st superdiagonal
   def d1(n):
      return -1.0/(2*n)

   ds = [d_1, d1]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 1)


def i2(nx, bc, coeff = 1.0):
   """Create operator for 2nd integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-2,3,2)
   nzrow = 1

   # Generate 2nd subdiagonal
   def d_2(n):
      return 1.0/(4*n*(n - 1))

   # Generate main diagonal
   def d0(n):
      return -1.0/(2*(n - 1)*(n + 1))

   # Generate 2nd superdiagonal
   def d2(n):
      return 1.0/(4*n*(n + 1))

   ds = [d_2, d0, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 2)


def i2lapl(nx, k, l, bc, coeff = 1.0):
   """Create operator for 2nd integral in x of Laplacian T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-2,3,2)
   nzrow = 1

   # Generate 2nd subdiagonal
   def d_2(n):
      return -(k**2 + l**2)/(4*n*(n - 1))

   # Generate main diagonal
   def d0(n):
      return (k**2 + l**2 + 2*n**2 - 2)/(2*(n - 1)*(n + 1))

   # Generate 2nd superdiagonal
   def d2(n):
      return -(k**2 + l**2)/(4*n*(n + 1))

   ds = [d_2, d0, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 2)


def i2laplh(nx, k, bc, coeff = 1.0):
   """Create operator for 2nd integral in x of horizontal Laplacian T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-2,3,2)
   nzrow = 1

   # Generate 2nd subdiagonal
   def d_2(n):
      return -k**2/(4*n*(n - 1))

   # Generate main diagonal
   def d0(n):
      return (k**2 + 2*n**2 - 2)/(2*(n - 1)*(n + 1))

   # Generate 2nd superdiagonal
   def d2(n):
      return -k**2/(4*n*(n + 1))

   ds = [d_2, d0, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 2)


def i4(nx, bc, coeff = 1.0):
   """Create operator for 4th integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3
   
   # Generate 4th subdiagonal
   def d_4(n):
      return 1.0/(16*n*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return -1.0/(4*n*(n - 3)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return 3.0/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return -1.0/(4*n*(n - 1)*(n + 1)*(n + 3)) 

   # Generate 4th superdiagonal
   def d4(n):
      return 1.0/(16*n*(n + 1)*(n + 2)*(n + 3))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def i4d2(nx, bc, coeff = 1.0):
   """Create operator for 4th integral in x of D_x^2 T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-2,3,2)
   nzrow = 3

   # Generate 2nd subdiagonal
   def d_2(n):
      return 1.0/(4*n*(n - 1))

   # Generate main diagonal
   def d0(n):
      return -1.0/(2*(n - 1)*(n + 1))

   # Generate 2nd superdiagonal
   def d2(n):
      return 1.0/(4*n*(n + 1))

   ds = [d_2, d0, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def i4d4(nx, bc, coeff = 1.0):
   """Create a quasi identity block of order 4"""

   return qid(nx,4, bc, coeff)


def i4lapl(nx, k, l, bc, coeff = 1.0):
   """Create operator for 4th integral in x of Laplacian T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3
   
   # Generate 4th subdiagonal
   def d_4(n):
      return -(k**2 + l**2)/(16*n*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return (k**2 + l**2 + n**2 - 2*n - 3)/(4*n*(n - 3)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return -(3*k**2 + 3*l**2 + 4*n**2 - 16)/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return (k**2 + l**2 + n**2 + 2*n - 3)/(4*n*(n - 1)*(n + 1)*(n + 3))

   # Generate 4th superdiagonal
   def d4(n):
      return -(k**2 + l**2)/(16*n*(n + 1)*(n + 2)*(n + 3))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def i4laplh(nx, k, bc, coeff = 1.0):
   """Create operator for 4th integral in x of horizontal Laplacian T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3
   
   # Generate 4th subdiagonal
   def d_4(n):
      return -k**2/(16*n*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return (k**2 + n**2 - 2*n - 3)/(4*n*(n - 3)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return -(3*k**2 + 4*n**2 - 16)/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return (k**2 + n**2 + 2*n - 3)/(4*n*(n - 1)*(n + 1)*(n + 3))

   # Generate 4th superdiagonal
   def d4(n):
      return -k**2/(16*n*(n + 1)*(n + 2)*(n + 3))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def i4lapl2(nx, k, l, bc, coeff = 1.0):
   """Create operator for 4th integral in x of Laplacian^2 T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3

   # Generate 4th subdiagonal
   def d_4(n):
      return (k**2 + l**2)**2/(16*n*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return -(k**2 + l**2)*(k**2 + l**2 + 2*n**2 - 4*n - 6)/(4*n*(n - 3)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return (3*k**4 + 6*k**2*l**2 + 8*k**2*n**2 - 32*k**2 + 3*l**4 + 8*l**2*n**2 - 32*l**2 + 8*n**4 - 40*n**2 + 32)/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return -(k**2 + l**2)*(k**2 + l**2 + 2*n**2 + 4*n - 6)/(4*n*(n - 1)*(n + 1)*(n + 3)) 

   # Generate 4th superdiagonal
   def d4(n):
      return (k**2 + l**2)**2/(16*n*(n + 1)*(n + 2)*(n + 3))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def i4lapl2h(nx, k, bc, coeff = 1.0):
   """Create operator for 4th integral in x of horizontal Laplacian^2 T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3

   # Generate 4th subdiagonal
   def d_4(n):
      return k**4/(16*n*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return -k**2*(k**2 + 2*n**2 - 4*n - 6)/(4*n*(n - 3)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return (3*k**4 + 8*k**2*n**2 - 32*k**2 + 8*n**4 - 40*n**2 + 32)/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return -k**2*(k**2 + 2*n**2 + 4*n - 6)/(4*n*(n - 1)*(n + 1)*(n + 3)) 

   # Generate 4th superdiagonal
   def d4(n):
      return k**4/(16*n*(n + 1)*(n + 2)*(n + 3))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, 4)


def qid(nx, q, bc, coeff = 1.0):
   """Create a quasi identity block of order q"""

   offsets = [0]
   diags = [[0]*q + [1]*(nx-q)]

   mat = coeff*spsp.diags(diags, offsets)
   return c1dbc.constrain(mat, bc, q)
