"""Module provides functions to generate sparse operators in a cartesian box with a single periodic dimension."""

from __future__ import division

import scipy.sparse as spsp
import cartesian_1d as c1d
import cartesian_boundary_2d as c2dbc

def convert_bc(bc):
   """Convert boundary dictionary into x and z kronecker product boundaries"""

   if bc['x'][0] < 0:
      bcx = bc['x']
   else:
      bcx = [0]

   if bc['z'][0] < 0:
      bcz = bc['z']
   else:
      bcz = [0]

   return (bcx, bcz)


def zblk(nx, nz, qx, qz, bc):
   """Create a block of zeros"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.zblk(nz,qz,bcz),c1d.zblk(nx,qx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, qx, qz)


def i2j2d2d2(nx, nz, bc):
   """Create a quasi identity block of order 2,2"""

   return qid(nx,nz,2,2, bc)


def i2j2(nx, nz, bc):
   """Create operator for 2nd integral in x,z of T_n(x)T_n(z)"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.i2(nz,bcz), c1d.i2(nx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, 2, 2)


def i2j2lapl(nx, nz, k, bc):
   """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.i2(nz,bcz), c1d.i2laplh(nx,k,bcz)) + spsp.kron(c1d.i2d2(nz,bcz), c1d.i2(nx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, 2, 2)


def i4j4(nx, nz, bc):
   """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4(nx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, 4, 4)


def i4j4lapl(nx, nz, k, bc):
   """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4laplh(nx,k,bcx)) + spsp.kron(c1d.i4d2(nz,bcz), c1d.i4(nx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, 4, 4)


def i4j4lapl2(nx, nz, k, bc):
   """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4lapl2h(nx,k,bcx)) + 2*spsp.kron(c1d.i4d2(nz,bcz), c1d.i4laplh(nx,k,bcx)) + spsp.kron(c1d.i4d4(nz,bcz), c1d.i4(nx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, 4, 4)


def qid(nx, nz, qx, qz, bc):
   """Create a quasi identity block order qx in x in qz in z"""

   bcx, bcz = convert_bc(bc)
   mat = spsp.kron(c1d.qid(nz,qz,bcz), c1d.qid(nx,qx,bcx))
   return c2dbc.constrain(mat, nx, nz, bc, qx, qz)
