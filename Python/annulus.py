"""Module provides functions to generate sparse operators in a cylindrical annulus."""

from __future__ import division

import scipy.sparse as spsp


def zblk(nr,nz):
   """Create a block of zeros"""

   return spsp.coo_matrix((nr*nz,nr*nz))


def i2j2x2(nr, nz, m, a, b):
   """Create operator for 2nd integral of x^2 T_n(x).""" 

   return spsp.identity(nr*nz)


def i2j2x2lapl(nr, nz, m, a, b):
   """Create operator for 2nd integral of x^2 Laplacian T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4(nr, nz, m, a, b):
   """Create operator for 4th integral of x^4 T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4lapl(nr, nz, m, a, b):
   """Create operator for 4th integral of x^4 Laplacian T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4lapl2(nr, nz, m, a, b):
   """Create operator for 4th integral of x^4 Laplacian^2 T_n(x).""" 

   return spsp.identity(nr*nz)
