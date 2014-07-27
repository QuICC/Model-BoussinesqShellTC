"""Module provides functions to transform chebyshev expansions for the radius of a shell between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import numpy.polynomial.legendre as leg
import scipy.special as spe

def rgrid(nr, a, b):
    """Create the radial Chebyshev grid"""

    return a*np.cos(np.pi*(np.arange(0,nr)+0.5)/nr) + b

def tgrid(maxl, m):
    """Create the theta Gauss-Legendre grid"""

    nt = 3*(maxl - m + 1)//2
    x, tmp = leg.leggauss(nt)

    return x

def torphys(spec):
    """Transform radial spectral coefficients to physical values"""

    n = len(spec)

    return fftpack.dct(spec,3)

def torcheb(phys):
    """Transform radial physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def totphys(spec, maxl, m):
    """Tranform theta spectral coefficients to physical values"""

    x = tgrid(maxl, m)
    mat = plm(maxl, m)
    phys = mat.dot(spec)

    return phys

def totleg(phys, maxl, m):
    """Tranform theta physical values to spectral coefficients"""

    nt = 3*(maxl - m + 1)//2
    x, w = leg.leggauss(nt)
    mat = plm(maxl, m).T
    spec = mat.dot(np.diag(w).dot(phys))

    return spec

def torphys2D(spec):
    """Transform 2D R spectral coefficients to 2D R physical values"""

    phys = np.zeros(spec.shape)
    for j in range(spec.shape[1]):
        phys[:,j] = torphys(spec[:,j])
    
    return phys

def torcheb2D(phys):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = torcheb(phys[:,j])

    return phys

def plm(maxl, m):
    """Compute the normalized associated legendre polynomial projection matrix"""

    x = tgrid(maxl, m)
    mat = np.zeros((len(x), maxl - m + 1))
    mat[:,0] = pmm(maxl, m)[:,0]
    mat[:,1] = np.sqrt(2.0*m + 3.0)*x*mat[:,0]
    for i, l in enumerate(range(m+2, maxl + 1)):
        mat[:,i+2] = np.sqrt((2.0*l + 1)/(l - m))*np.sqrt((2.0*l - 1.0)/(l + m))*x*mat[:,i+1] - np.sqrt((2.0*l + 1)/(2.0*l - 3.0))*np.sqrt((l + m - 1.0)/(l + m))*np.sqrt((l - m - 1.0)/(l - m))*mat[:,i]

    return mat

def pmm(maxl, m):
    """Compute the normalized associated legendre polynomial of order and degree m"""

    x = tgrid(maxl, m)
    mat = np.zeros((len(x), 1))
    mat[:,0] = 1.0/np.sqrt(2.0)
    sx = np.sqrt(1.0 - x**2)
    for i in range(1, m+1):
        mat[:,0] = -np.sqrt((2.0*i + 1.0)/(2.0*i))*sx*mat[:,0]

    return mat
