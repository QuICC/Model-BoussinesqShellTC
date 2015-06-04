"""Module provides functions to transform chebyshev expansions for the radius of a shell between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np
import numpy.polynomial.legendre as leg
import scipy.special as spe

min_r_points = 1000
min_th_points = 1000
min_phi_points = 1000

def rgrid(nr, a, b):
    """Create the radial Chebyshev grid"""

    gN = max(min_r_points, nr)

    return a*np.cos(np.pi*(np.arange(0,gN)+0.5)/gN) + b

def tgrid(maxl, m):
    """Create the theta Gauss-Legendre grid"""

    nt = max(min_th_points,3*(maxl - m + 1)//2)
    nt = nt + (nt+1)%2
    x, tmp = leg.leggauss(nt)

    return x

def phgrid(m):
    """Create a phi grid for given harmonic order"""

    return np.linspace(0, 2*np.pi, max(min_phi_points,3*m))

def grid_2d(nr, a, b, maxl, m):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr, a, b)
    th = np.arccos(tgrid(maxl, m))
    rmesh, thmesh = np.meshgrid(r, th)
    X = rmesh * np.sin(thmesh)
    Y = rmesh * np.cos(thmesh)

    return (X, Y)

def grid_eq(nr, a, b, m):
    """Compute the 2D grid for the equatorial contours"""

    r = rgrid(nr, a, b)
    phi = phgrid(m)
    rmesh, phimesh = np.meshgrid(r, phi)
    X = rmesh * np.cos(phimesh)
    Y = rmesh * np.sin(phimesh)

    return (X, Y)

def torphys(spec):
    """Transform radial spectral coefficients to physical values"""

    if len(spec) < min_r_points:
        data = np.hstack((spec, np.zeros(min_r_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data,3)

def torcheb(phys):
    """Transform radial physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

def totphys(spec, maxl, m):
    """Tranform theta spectral coefficients to physical values"""
    
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

    phys = np.zeros((max(spec.shape[0],min_r_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:,j] = torphys(spec[:,j])
    
    return phys

def torcheb2D(phys):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:,j] = torcheb(phys[:,j])

    return phys

def toslice(spec, nr, maxl, m):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, maxl-m+1), order = 'F')

    rphys = torphys2D(spec)

    rspec = np.transpose(rphys)
    phys = totphys(rspec, maxl, m)
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
