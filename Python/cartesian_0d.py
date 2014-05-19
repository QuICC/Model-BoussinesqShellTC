"""Module provides functions to generate sparse operators in a triply periodic cartesian box."""

from __future__ import division

import scipy.sparse as spsp


def id(bc):
   """Create an identity block"""

   return spsp.identity(1)


def lapl(k, l, m, bc):
   """Create operator for triply periodic Laplacian"""

   return spsp.identity(1)


def lapl2(k, l, m, bc):
   """Create operator for triply periodic Laplacian^2"""

   return spsp.identity(1)
