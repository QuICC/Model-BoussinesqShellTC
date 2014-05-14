"""Module provides functions to generate sparse operators in a cartesian box with two periodic directions."""

import scipy.sparse as spsp
import utils


def zblk(nx):
   """Create a block of zeros"""

   return spsp.coo_matrix((nx,nx))


def id(nx):
   """Create an identity block"""

   return spsp.identity(nx)


def i1(nx):
   """Create operator for 1st integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-1,2,2)
   nzrow = 0

   # Generate 1st subdiagonal
   def d_1(n):
      return -1.0/(2*(n-1))

   # Generate 1st superdiagonal
   def d1(n):
      return 1.0/(2*(n+1))

   ds = [d_1, d1]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i2(nx):
   """Create operator for 2nd integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-2,3,2)
   nzrow = 1

   # Generate 2nd subdiagonal
   def d_2(n):
      return 1/(4*(n - 2)*(n - 1))

   # Generate main diagonal
   def d0(n):
      return -1/(2*(n - 1)*(n + 1))

   # Generate 2nd superdiagonal
   def d2(n):
      return 1/(4*(n + 1)*(n + 2))

   ds = [d_2, d0, d2]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i2lapl(nx, k, l):
   """Create operator for 2nd integral in x of Laplacian T_n(x)"""

   return spsp.identity(nx)


def i4(nx):
   """Create operator for 4th integral in x of T_n(x)"""

   ns = np.arange(0, nx, 1)
   offsets = np.arange(-4,5,2)
   nzrow = 3
   
   # Generate 4th subdiagonal
   def d_4(n):
      return 1/(16*(n - 4)*(n - 3)*(n - 2)*(n - 1))

   # Generate 2nd subdiagonal
   def d_2(n):
      return -1/(4*(n - 3)*(n - 2)*(n - 1)*(n + 1)) 

   # Generate main diagonal
   def d0(n):
      return 3/(8*(n - 2)*(n - 1)*(n + 1)*(n + 2))

   # Generate 2nd superdiagonal
   def d2(n):
      return -1/(4*(n - 1)*(n + 1)*(n + 2)*(n + 3))

   # Generate 4th superdiagonal
   def d4(n):
      return 1/(16*(n + 1)*(n + 2)*(n + 3)*(n + 4))

   ds = [d_4, d_2, d0, d2, d4]
   diags = utils.build_diagonals(ns, nzrow, ds, offsets)

   return spsp.diags(diags, offsets)


def i4lapl(nx, k, l):
   """Create operator for 4th integral in x of Laplacian T_n(x)"""

   return spsp.identity(nx)


def i4lapl2(nx, k, l):
   """Create operator for 4th integral in x of Laplacian^2 T_n(x)"""

   return spsp.identity(nx)
