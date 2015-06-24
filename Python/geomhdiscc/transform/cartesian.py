"""Module provides functions to transform chebyshev expansions between physical and spectral space."""

from __future__ import division
from __future__ import unicode_literals

import scipy.fftpack as fftpack
import numpy as np

from geomhdiscc.transform.spherical import eqgrid

def grid(nx):
    """Create the Chebyshev grid"""

    return np.cos(np.pi*(np.arange(0,nx)+0.5)/nx)

grid_fast = grid
grid_slow = grid

def grid_2d(nx, nz):
    """Create 2D grid"""

    xGrid = grid(nx)
    zGrid = grid(nz)
    xMesh, yMesh = np.meshgrid(xGrid, zGrid)

    return(xMesh, yMesh)

def grid_fast_per(nx, m):
    """Create 2D grid"""

    phi = eqgrid(m)
    xGrid = grid(nx)
    xMesh, yMesh = np.meshgrid(xGrid, phi)

    return(xMesh, yMesh)

def grid_slow_per(nz, m):
    """Create 2D grid"""

    phi = eqgrid(m)
    zGrid = grid(nz)
    xMesh, yMesh = np.meshgrid(zGrid, phi)

    return(xMesh, yMesh)

def tophys(spec):
    """Transform R spectral coefficients to physical values"""

    n = len(spec)
    phys = np.zeros(spec.shape, dtype = spec.dtype)
    if spec.dtype == 'complex_':
        phys.real = fftpack.dct(spec.real,3)
        phys.imag = fftpack.dct(spec.imag,3)
    else:
        phys = fftpack.dct(spec,3)

    return phys

def tocheb(phys):
    """Transform physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys,2)/(2*n)

grid_1d = grid
toprofile = tophys

def tophys2d(spec):
    """Transform 2D spectral coefficients to 2D physical values"""

    phys = spec.copy()

    for i in range(spec.shape[0]):
        phys[i,:] = tophys(spec[i,:])
    for j in range(spec.shape[1]):
        phys[:,j] = tophys(phys[:,j])

    phys = phys.T
    
    return phys

def toslice(spec, nx, nz):
    """Convert linear spectral to physical values on grid"""

    spec = np.reshape(spec, (nx, nz), order ='F')

    phys = tophys2d(spec)
    return phys

def tocheb2d(phys):
    """Transform 2D physical array to 2D spectral coefficients"""

    spec = phys.copy()

    for j in range(phys.shape[1]):
        spec[:,j] = tocheb(phys[:,j])
    for i in range(phys.shape[0]):
        spec[i,:] = tocheb(spec[i,:])
    
    return spec

def tophys3d(spec):
    """Transform 3D spectral coefficients to 3D physical values"""

    phys = spec.copy()

    for i in range(spec.shape[1]):
        for j in range(spec.shape[2]):
            phys[:,i,j] = tophys(phys[:,i,j])
    for i in range(spec.shape[0]):
        for j in range(spec.shape[2]):
            phys[i,:,j] = tophys(phys[i,:,j])
    for i in range(spec.shape[0]):
        for j in range(spec.shape[1]):
            phys[i,j,:] = tophys(spec[i,j,:])
    
    return phys

def tovolume(spec, nx, ny, nz):
    """Convert linear spectral to physical values on grid"""

    spec = np.reshape(spec, (nx, ny, nz), order ='F')

    phys = tophys3d(spec)
    return phys

def tocheb3d(phys):
    """Transform 3D physical array to 3D spectral coefficients"""

    for i in range(phys.shape[1]):
        for j in range(phys.shape[2]):
            phys[:,i,j] = tocheb(phys[:,i,j])
    for i in range(phys.shape[0]):
        for j in range(phys.shape[2]):
            phys[i,:,j] = tocheb(phys[i,:,j])
    for i in range(phys.shape[0]):
        for j in range(phys.shape[1]):
            phys[i,j,:] = tocheb(phys[i,j,:])
    
    return phys
