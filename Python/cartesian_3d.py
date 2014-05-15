"""Module provides functions to generate sparse operators in a cartesian box."""

import scipy.sparse as spsp
import cartesian_1d as c1d
import cartesian_2d as c2d


def zblk(nx,ny,nz):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx*nx*nz,nx*ny*nz))


def i2j2k2d2d2d2(nx, ny, nz):
   """Create a quasi identity block of order 2,2,2"""

   return qid(nx,ny,nz,2,2,2)


def i2j2k2(nx, ny, nz):
   """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

   return spsp.kron(c1d.i2(ny), c2d.i2j2(nx,nz))


def i2j2k2lapl(nx, ny, nz):
   """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

   op = spsp.kron(c1d.i2(ny),spsp.kron(c1d.i2(nz),c1d.i2d2(nx)))
   op = op + spsp.kron(c1d.i2d2(ny),spsp.kron(c1d.i2(nz),c1d.i2(nx)))
   op = op + spsp.kron(c1d.i2(ny),spsp.kron(c1d.i2d2(nz),c1d.i2(nx)))

   return  op


def i4j4k4(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

   return spsp.kron(c1d.i4(ny), c2d.i4j4(nx,nz))


def i4j4k4lapl(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

   op = spsp.kron(c1d.i4(ny),spsp.kron(c1d.i4(nz),c1d.i4d2(nx)))
   op = op + spsp.kron(c1d.i4d2(ny),spsp.kron(c1d.i4(nz),c1d.i4(nx)))
   op = op + spsp.kron(c1d.i4(ny),spsp.kron(c1d.i4d2(nz),c1d.i4(nx)))

   return  op


def i4j4k4lapl2(nx, ny, nz):
   """Create operator for 4th integral in x,y,z of Laplacian^2 T_n(x)T_n(y)T_n(z)"""

   op = spsp.kron(c1d.i4(ny),spsp.kron(c1d.i4(nz),c1d.i4d4(nx)))
   op = op + spsp.kron(c1d.i4d4(ny),spsp.kron(c1d.i4(nz),c1d.i4(nx)))
   op = op + spsp.kron(c1d.i4(ny),spsp.kron(c1d.i4d4(nz),c1d.i4(nx)))
   op = op + 2*spsp.kron(c1d.i4d2(ny),spsp.kron(c1d.i4(nz),c1d.i4d2(nx)))
   op = op + 2*spsp.kron(c1d.i4d2(ny),spsp.kron(c1d.i4d2(nz),c1d.i4(nx)))
   op = op + 2*spsp.kron(c1d.i4(ny),spsp.kron(c1d.i4d2(nz),c1d.i4d2(nx)))

   return  op


def qid(nx, ny, nz, qx, qy, qz):
   """Create a quasi identity block order qx,qy,qz"""

   return spsp.kron(c1d.qid(ny,qy), c2d.qid(nx,nz,qx,qz))
