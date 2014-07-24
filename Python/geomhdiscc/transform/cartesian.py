"""Module provides functions to transform chebyshev expansions between physical and spectral space."""

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

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    for i in range(phys.shape[0]):
        phys[i,:] = tocheb(phys[i,:])
    for j in range(phys.shape[1]):
        phys[:,j] = tocheb(phys[:,j])
    
    return phys
