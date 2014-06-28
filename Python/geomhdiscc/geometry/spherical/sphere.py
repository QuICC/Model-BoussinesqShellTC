"""Module provides functions to generate sparse operators in a sphere."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.spherical_boundary as sphbc


def zblk(nr):
   """Create a block of zeros"""

   return spsp.lil_matrix((nr,nr))


def i2x2(nr, l, m):
   """Create operator for 2nd integral of x^2 T_n(x)."""

   parity = l%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-2,3)
   nzrow = 1

   # Generate 2nd subdiagonal
   def d_2(n):
      return 1.0/(16.0*n*(n - 1.0))

   # Generate 1st subdiagonal
   def d_1(n):
      return -1.0/(8.0*n*(n - 1.0)*(n + 1.0))

   # Generate main diagonal
   def d0(n):
      return -1.0/(8.0*(n - 1.0)*(n + 1.0))

   # Generate 1st superdiagonal
   def d1(n):
      return -d1(n)

   # Generate 2nd superdiagonal
   def d2(n):
      return d_2(n + 1.0)

   ds = [d_2, d_1, d0, d1, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i2x2lapl(nr, l, m):
   """Create operator for 2nd integral of x^2 Laplacian T_n(x).""" 

   parity = l%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-1,2)
   nzrow = 1
   
   # Generate 1st subdiagonal
   def d_1(n):
      return  -((l - n + 2.0)*(l + n - 1.0))/(4.0*n*(n - 1.0))

   # Generate main diagonal
   def d0(n):
      return (l**2 + l + n**2 - 1.0)/(2.0*(n - 1.0)*(n + 1.0))

   # Generate 1st superdiagonal
   def d1(n):
      return d_1(n + 1.0) + (2*n + 1.0)/(2.0*n*(n + 1.0))

   ds = [d_1, d0, d1]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i4x4(nr, l, m):
   """Create operator for 4th integral of x^4 T_n(x).""" 
   
   parity = l%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-4,5)
   nzrow = 3
   
   # Generate 4th subdiagonal
   def d_4(n):
      return 1.0/(256.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

   # Generate 3rd subdiagonal
   def d_3(n):
      return 3.0/(64.0*n*(n - 1.0)*(n + 1.0)*(n - 2.0)*(n - 3.0)) 

   # Generate 2nd subdiagonal
   def d_2(n):
      return -(n**2 - 19.0)/(64.0*n*(n - 1.0)*(n + 1.0)*(n - 2.0)*(n + 2.0)*(n - 3.0)) 

   # Generate 1st subdiagonal
   def d_1(n):
      return -(3.0*(3.0*n + 11.0))/(64.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n - 3.0)*(n + 3.0)) 

   # Generate main diagonal
   def d0(n):
      return (3.0*(n**2 - 29.0))/(128.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

   # Generate 1st superdiagonal
   def d1(n):
      return (3.0*(3.0*n - 11.0))/(64.0*n*(n - 1.0)*(n + 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)) 

   # Generate 2nd superdiagonal
   def d2(n):
      return -(n**2 - 19.0)/(64.0*n*(n - 1.0)*(n + 1.0)*(n - 2.0)*(n + 2.0)*(n + 3.0))

   # Generate 3rd superdiagonal
   def d3(n):
      return -d_3(n + 2.0)

   # Generate 4th superdiagonal
   def d4(n):
      return d_4(n + 3.0)

   ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i4x4lapl(nr, l, m):
   """Create operator for 4th integral of x^4 Laplacian T_n(x).""" 

   parity = l%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-3,4)
   nzrow = 3
   
   # Generate 3rd subdiagonal
   def d_3(n):
      return -((l - n + 6.0)*(l + n - 5.0))/(64.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

   # Generate 2nd subdiagonal
   def d_2(n):
      return (l**2*n - 5.0*l**2 + l*n - 5.0*l + n**3 - 3.0*n**2 - 28.0*n + 96.0)/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0)) 

   # Generate 1st subdiagonal
   def d_1(n):
      return (l**2*n + 17.0*l**2 + l*n + 17.0*l - n**3 + 24.0*n**2 - 5.0*n - 294.0)/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 2.0)*(n + 1.0)) 

   # Generate main diagonal
   def d0(n):
      return -(l**2*n**2 - 19.0*l**2 + l*n**2 - 19.0*l + n**4 - 37.0*n**2 + 312.0)/(16.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0)) 

   # Generate 1st superdiagonal
   def d1(n):
      return (l**2*n - 17.0*l**2 + l*n - 17.0*l - n**3 - 24.0*n**2 - 5.0*n + 294.0)/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 1.0)) 

   # Generate 2nd superdiagonal
   def d2(n):
      return (l**2*n + 5.0*l**2 + l*n + 5.0*l + n**3 + 3.0*n**2 - 28.0*n - 96.0)/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0)) 

   # Generate 3rd superdiagonal
   def d3(n):
      return -((l - n - 5.0)*(l + n + 6.0))/(64.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0)) 

   ds = [d_3, d_2, d_1, d0, d1, d2, d3]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i4x4lapl2(nr, l, m):
   """Create operator for 4th integral of x^4 Laplacian^2 T_n(x).""" 
   
   parity = l%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-2,3)
   nzrow = 3
   
   # Generate 2nd subdiagonal
   def d_2(n):
      return ((l + n - 5.0)*(l - n + 4.0)*(l - n + 6.0)*(l + n - 3.0))/(16.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

   # Generate 1st subdiagonal
   def d_1(n):
      return -((l + n - 3.0)*(l - n + 4.0)*(l**2 + l + n**2 - 2.0*n - 9.0))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

   # Generate main diagonal
   def d0(n):
      return (3.0*l**4 + 6.0*l**3 + 2.0*l**2*n**2 - 47.0*l**2 + 2.0*l*n**2 - 50.0*l + 3.0*n**4 - 51.0*n**2 + 228.0)/(8.0*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0)) 

   # Generate 1st superdiagonal
   def d1(n):
      return -((l + n + 4.0)*(l - n - 3.0)*(l**2 + l + n**2 + 2.0*n - 9.0))/(4.0*n*(n + 3.0)*(n - 1.0)*(n + 1.0)) 

   # Generate 2nd superdiagonal
   def d2(n):
      return ((l - n - 5.0)*(l - n - 3.0)*(l + n + 6.0)*(l + n + 4.0))/(16.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0)) 

   ds = [d_2, d_1, d0, d1, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def block_diag(nr, nl, m, op):
   """Create a block diagonal matrix of operator."""

   B = op(nr,m,m)
   for l in np.arange(m+1,nl-1):
      A = op(nr,l,m)
      B = spsp.block_diag((B,A))

   return B


