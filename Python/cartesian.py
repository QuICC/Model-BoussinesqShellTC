"""Module provides functions to generate sparse operators in a cartesian box."""

import scipy.sparse as spsp


##################################################
# 3D Cartesian box

def zblk(nx,ny,nz):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx*nx*nz,nx*ny*nz))


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

##################################################
# 2D Cartesian box with single periodic direction

def zblk(nx,ny):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx*nx,nx*ny))


def i2j2(nx, ny, k):
   """Create operator for 2nd integral in x,z of T_n(x)T_n(z)"""

   return spsp.identity(nx*ny)


def i2j2lapl(nx, ny, k):
   """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

   return spsp.identity(nx*ny)


def i4j4(nx, ny, k):
   """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

   return spsp.identity(nx*ny)


def i4j4lapl(nx, ny, k):
   """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

   return spsp.identity(nx*ny)


def i4j4lapl2(nx, ny, k):
   """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

   return spsp.identity(nx*ny)


##################################################
# 1D Cartesian box with two periodic directions

def zblk(nx):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx,nx))


def i2(nx, k, l):
   """Create operator for 2nd integral in x of T_n(x)"""

   return spsp.identity(nx)


def i2lapl(nx, k, l):
   """Create operator for 2nd integral in x of Laplacian T_n(x)"""

   return spsp.identity(nx)


def i4(nx, k, l):
   """Create operator for 4th integral in x of T_n(x)"""

   return spsp.identity(nx)


def i4lapl(nx, k, l):
   """Create operator for 4th integral in x of Laplacian T_n(x)"""

   return spsp.identity(nx)


def i4lapl2(nx, k, l):
   """Create operator for 4th integral in x of Laplacian^2 T_n(x)"""

   return spsp.identity(nx)
