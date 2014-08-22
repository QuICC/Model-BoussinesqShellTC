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

def torphys(spec):
    """Transform R spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def torcheb(phys):
    """Transform R physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tozphys(spec):
    """Transform Z spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def tozcheb(phys):
    """Transform Z physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def tophys2d(spec):
    """Transform 2D spectral coefficients to 2D physical values"""

    phys = spec.copy()

    for i in range(spec.shape[0]):
        phys.real[i,:] = tozphys(spec.real[i,:])
        phys.imag[i,:] = tozphys(spec.imag[i,:])
    for j in range(spec.shape[1]):
        phys.real[:,j] = torphys(phys.real[:,j])
        phys.imag[:,j] = torphys(phys.imag[:,j])
    
    return phys

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    spec = phys.copy()
    for j in range(phys.shape[1]):
        spec[:,j] = torcheb(phys[:,j])
    for i in range(phys.shape[0]):
        spec[i,:] = tozcheb(spec[i,:])
    
    return spec
