"""Module provides functions to work with chebyshev expansions. It is mainly meant for testing purposes."""

import scipy as sp
import scipy.fftpack as ff
import numpy as np

def grid(nx):
   """Create the Chebyshev grid"""

   return np.cos(np.pi*np.linspace(0,1,nx))

def tophys(spec):
   """Transform spectral coefficients to physical values"""

   n = len(spec)

   if n > 2: 
      t = spec[1:-1]
      ifft = np.fft.ifft((n-1)*np.append(spec, t[::-1]))

      return ifft.real[:n]
   else:
      return spec/2

def tocheb(phys):
   """Transform physical values to spectral coefficients"""

   n = len(phys)

   if n > 2: 
      t = phys[1:-1]
      fft = np.fft.fft(np.append(phys, t[::-1]))

      return fft.real[:n]/(n-1)
   else:
      return phys*2

