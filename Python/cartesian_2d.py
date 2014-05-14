"""Module provides functions to generate sparse operators in a cartesian box with a single periodic dimension."""

import scipy.sparse as spsp
import cartesian_1d as c1d


def zblk(nx, nz):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx*nx,nx*nz))


def i2j2d2d2(nx, nz):
   """Create a quasi identity block"""

   return spsp.kron(c1d.i2d2(nz), c1d.i2d2(nx))


def i2j2(nx, nz):
   """Create operator for 2nd integral in x,z of T_n(x)T_n(z)"""

   return spsp.kron(c1d.i2(nz), c1d.i2(nx))


def i2j2lapl(nx, nz, k):
   """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

   return spsp.kron(c1d.i2(nz), c1d.i2laplh(nx,k)) + spsp.kron(c1d.i2d2(nz), c1d,i2(nx))


def i4j4(nx, nz):
   """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

   return spsp.kron(c1d.i4(nz), c1d.i4(nx))


def i4j4lapl(nx, nz, k):
   """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

   return spsp.identity(nx*nz)


def i4j4lapl2(nx, nz, k):
   """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

   return spsp.identity(nx*nz)
