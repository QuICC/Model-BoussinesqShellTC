"""Module provides functions to transform chebyshev expansions for the radius of an annulus between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np


def rgrid(nr, a, b):
    """Create the radial Chebyshev grid"""

    return a*np.cos(np.pi*(np.arange(0,nr)+0.5)/nr) + b

def zgrid(nz):
    """Create the z Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,nz)+0.5)/nz)

def tophys(spec):
    """Transform spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def tocheb(phys):
    """Transform physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tophys2d(spec):
    """Transform 2D spectral coefficients to 2D physical values"""

    for i in range(spec.shape[0]):
        spec.real[i,:] = tophys(spec.real[i,:])
        spec.imag[i,:] = tophys(spec.imag[i,:])
    for j in range(spec.shape[1]):
        spec.real[:,j] = tophys(spec.real[:,j])
        spec.imag[:,j] = tophys(spec.imag[:,j])
    
    return spec

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    for i in range(phys.shape[0]):
        phys[i,:] = tocheb(phys[i,:])
    for j in range(phys.shape[1]):
        phys[:,j] = tocheb(phys[:,j])
    
    return phys
