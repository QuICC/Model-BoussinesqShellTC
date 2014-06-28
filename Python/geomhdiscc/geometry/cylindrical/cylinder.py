"""Module provides functions to generate sparse operators in a cylinder."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cylindrical.cylinder_radius as rad
import geomhdiscc.geometry.cylindrical.cylindrical_boundary as cylbc


def convert_bc(bc):
   """Convert boundary dictionary into r and z kronecker product boundaries"""

   if bc['r'][0] < 0:
      bcr = bc['r']
   else:
      bcr = [0]

   if bc['z'][0] < 0:
      bcz = bc['z']
   else:
      bcz = [0]

   return (bcr, bcz)


def zblk(nr,nz):
   """Create a block of zeros"""

   bcr, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.zblk(nz,qz,bcz),rad.zblk(nr,qr,bcr))
   return cylbc.constrain(mat, nr, nz, bc, qr, qz)


def i2j2x2(nr, nz, m):
   """Create operator for 2nd integral of x^2 T_n(x).""" 

   return spsp.identity(nr*nz)


def i2j2x2lapl(nr, nz, m):
   """Create operator for 2nd integral of x^2 Laplacian T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4(nr, nz, m):
   """Create operator for 4th integral of x^4 T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4lapl(nr, nz, m):
   """Create operator for 4th integral of x^4 Laplacian T_n(x).""" 

   return spsp.identity(nr*nz)


def i4j4x4lapl2(nr, nz, m):
   """Create operator for 4th integral of x^4 Laplacian^2 T_n(x).""" 

   return spsp.identity(nr*nz)

