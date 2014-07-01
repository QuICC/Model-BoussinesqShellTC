"""Module provides functions to generate sparse operators for the radial direction in a cylinder."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.cylindrical_radius_boundary as cylbc


def zblk(nr, q, bc):
   """Create a block of zeros"""

   mat = spsp.lil_matrix((nr,nr))
   return cylbc.constrain(mat,bc,q)


def i2x2(nr, m, bc, coeff = 1.0):
   """Create operator for 2nd integral of x^2 T_n(x)."""

   parity = m%2
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

   mat = coeff*spsp.diags(diags, offsets)
   return cylbc.constrain(mat, bc, 1)


def i2x2lapl(nr, m, bc, coeff = 1.0):
   """Create operator for 2nd integral of x^2 Laplacian T_n(x).""" 

   parity = m%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-1,2)
   nzrow = 1
   
   # Generate 1st subdiagonal
   def d_1(n):
      return  -(m - n + 2.0)*(m + n - 2.0)/(4.0*n*(n - 1.0))

   # Generate main diagonal
   def d0(n):
      return (m**2 + n**2 - 2.0)/(2.0*(n - 1.0)*(n + 1.0))

   # Generate 1st superdiagonal
   def d1(n):
      return -(m - n - 2.0)*(m + n + 2.0)/(4.0*n*(n + 1.0))

   ds = [d_1, d0, d1]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return cylbc.constrain(mat, bc, 1)


def i4x4(nr, m, bc, coeff = 1.0):
   """Create operator for 4th integral of x^4 T_n(x).""" 
   
   parity = m%2
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

   mat = coeff*spsp.diags(diags, offsets)
   return cylbc.constrain(mat, bc, 2)


def i4x4lapl(nr, m, bc, coeff = 1.0):
   """Create operator for 4th integral of x^4 Laplacian T_n(x).""" 

   parity = m%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-3,4)
   nzrow = 3
   
   # Generate 3rd subdiagonal
   def d_3(n):
      return -(m - n + 6.0)*(m + n - 6.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

   # Generate 2nd subdiagonal
   def d_2(n):
      return (n - 5.0)*(m**2 + n**2 + 2.0*n - 24.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)) 

   # Generate 1st subdiagonal
   def d_1(n):
      return (m**2*n + 17.0*m**2 - n**3 + 27.0*n**2 - 8.0*n - 372.0)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)) 

   # Generate main diagonal
   def d0(n):
      return -(m**2*n**2 - 19.0*m**2 + n**4 - 43.0*n**2 + 396.0)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0)) 

   # Generate 1st superdiagonal
   def d1(n):
      return (m**2*n - 17.0*m**2 - n**3 - 27.0*n**2 - 8.0*n + 372.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0)) 

   # Generate 2nd superdiagonal
   def d2(n):
      return (n + 5.0)*(m**2 + n**2 - 2.0*n - 24.0)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0)) 

   # Generate 3rd superdiagonal
   def d3(n):
      return -(m - n - 6.0)*(m + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0)) 

   ds = [d_3, d_2, d_1, d0, d1, d2, d3]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return cylbc.constrain(mat, bc, 2)


def i4x4lapl2(nr, m, bc, coeff = 1.0):
   """Create operator for 4th integral of x^4 Laplacian^2 T_n(x).""" 
   
   parity = m%2
   ns = np.arange(parity, 2*nr, 2)
   offsets = np.arange(-2,3)
   nzrow = 3
   
   # Generate 2nd subdiagonal
   def d_2(n):
      return (m - n + 4.0)*(m - n + 6.0)*(m + n - 6.0)*(m + n - 4.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

   # Generate 1st subdiagonal
   def d_1(n):
      return -(m - n + 4.0)*(m + n - 4.0)*(m**2 + n**2 - 2.0*n - 12.0)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

   # Generate main diagonal
   def d0(n):
      return (3.0*m**4 + 2.0*m**2*n**2 - 68.0*m**2 + 3.0*n**4 - 68.0*n**2 + 416.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)) 

   # Generate 1st superdiagonal
   def d1(n):
      return -(m - n - 4.0)*(m + n + 4.0)*(m**2 + n**2 + 2.0*n - 12.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0)) 

   # Generate 2nd superdiagonal
   def d2(n):
      return (m - n - 6.0)*(m - n - 4.0)*(m + n + 4.0)*(m + n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0)) 

   ds = [d_2, d_1, d0, d1, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   mat = coeff*spsp.diags(diags, offsets)
   return cylbc.constrain(mat, bc, 2)

