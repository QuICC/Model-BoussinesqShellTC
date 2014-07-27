"""Module provides functions to transform chebyshev expansions for the radius of a sphere between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np

def rgrid(nr):
    """Create the radial Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,2*nr)+0.5)/(2*nr))

def torphys(spec, parity):
    """Transform R spectral coefficients to R physical values"""

    n = 2*len(spec)
    full = np.array([0.0]*n)
    full[np.arange(parity,n,2)] = spec

    return fftpack.dct(full,3)

def torcheb(phys, parity):
    """Transform R physical values to R spectral coefficients"""

    n = len(phys)
    spec = fftpack.dct(phys,2)/(2*n)

    return spec[np.arange(parity,n,2)]
