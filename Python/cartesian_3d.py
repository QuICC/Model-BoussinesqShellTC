"""Module provides functions to generate sparse operators in a cartesian box."""

import scipy.sparse as spsp


def zblk(nx,ny,nz):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx*nx*nz,nx*ny*nz))


def id(nx, ny, nz):
   """Create an identity block"""

   return spsp.identity(nx*ny*nz)


def i2j2k2(nx, ny, nz):
   """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

   return spsp.identity(nx*ny*nz)


def i2j2k2lapl(nx, ny, nz):
   """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

   return spsp.identity(nx*ny*nz)


def i4j4k4(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

   return spsp.identity(nx*ny*nz)


def i4j4k4lapl(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

   return spsp.identity(nx*ny*nz)


def i4j4k4lapl2(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of Laplacian^2 T_n(x)T_n(y)T_n(z)"""

   return spsp.identity(nx*ny*nz)
