"""Module provides functions to transform chebyshev expansions between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np


def grid(nx):
    """Create the Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,nx)+0.5)/nx)

def tophys(spec):
    """Transform R spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def tocheb(phys):
    """Transform physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tophys2d(spec):
    """Transform 2D spectral coefficients to 2D physical values"""

    phys = spec.copy()

    if spec.dtype == 'complex_':
        for i in range(spec.shape[0]):
            phys.real[i,:] = tophys(spec.real[i,:])
            phys.imag[i,:] = tophys(spec.imag[i,:])
        for j in range(spec.shape[1]):
            phys.real[:,j] = tophys(phys.real[:,j])
            phys.imag[:,j] = tophys(phys.imag[:,j])
    else:
        for i in range(spec.shape[0]):
            phys[i,:] = tophys(spec[i,:])
        for j in range(spec.shape[1]):
            phys[:,j] = tophys(phys[:,j])
    
    return phys

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = tocheb(phys[:,j])
    for i in range(phys.shape[0]):
        phys[i,:] = tocheb(phys[i,:])
    
    return phys
