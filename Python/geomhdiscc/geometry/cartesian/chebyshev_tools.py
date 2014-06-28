"""Module provides functions to work with chebyshev expansions. It is mainly meant for testing purposes."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np

def grid(nx):
   """Create the Chebyshev grid"""

   return np.cos(np.pi*(np.arange(0,nx)+0.5)/nx)

def tophys(spec):
   """Transform spectral coefficients to physical values"""

   n = len(spec)

   return fftpack.dct(spec,3)

def tocheb(phys):
   """Transform physical values to spectral coefficients"""

   n = len(phys)

   return fftpack.dct(phys,2)/(2*n)

