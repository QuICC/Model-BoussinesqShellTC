"""Module provides functions to transform chebyshev expansions for the radius in a cylinder between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np

def grid(nr):
    """Create the Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,2*nr)+0.5)/(2*nr))

def tophys(spec, parity):
    """Transform spectral coefficients to physical values"""

    n = 2*len(spec)
    full = np.array([0.0]*n)
    full[np.arange(parity,n,2)] = spec

    return fftpack.dct(full,3)

def tocheb(phys, parity):
    """Transform physical values to spectral coefficients"""

    n = len(phys)
    spec = fftpack.dct(phys,2)/(2*n)

    return spec[np.arange(parity,n,2)]
