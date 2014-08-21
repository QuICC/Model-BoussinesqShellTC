"""Module provides functions to generate sparse operators in a sphere using a spherical harmonics expansion in the angular directions."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.geometry.spherical.sphere_radius as rad
import geomhdiscc.geometry.spherical.sphere_sh as sh
import geomhdiscc.geometry.spherical.sphere_boundary as sphbc


def convert_bc(bc):
    """Convert boundary condition to be suitable for kronecker product boundaries"""

    if bc[0] < 0:
        bcr = bc
    else:
        bcr = rad.radbc.no_bc()
        for key, val in bc['r'].items():
            if key != 0:
                bcr[key] = val

    return bcr

def sh_coeff(coeff):
    """Compute coefficients from spherical harmonic expansion"""

    if coeff == 'laplh':
        def fct(x):
            return x*(x + 1.0)
    else:
        def fct(x):
            return 1.0
    
    return fct

def zblk(nr, maxl, m, bc):
    """Create a block of zeros"""

    bcr = convert_bc(bc)

    mat = rad.zblk(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,rad.zblk(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i2x2(nr, maxl, m, bc, coeff = 1.0, with_sh_coeff = None):
    """Create a i2x2 radial operator kronecker with an identity"""
    
    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    mat = coeff*shc(m)*rad.i2x2(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x2(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i2x2lapl(nr, maxl, m, bc, coeff = 1.0, with_sh_coeff = None):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    mat = coeff*shc(m)*rad.i2x2lapl(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x2lapl(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i4x4(nr, maxl, m, bc, coeff = 1.0, with_sh_coeff = None):
    """Create a i4x4 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    mat = coeff*shc(m)*rad.i4x4(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i4x4lapl(nr, maxl, m, bc, coeff = 1.0, with_sh_coeff = None):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    bcr = convert_bc(bc)

    shc = sh_coeff(with_sh_coeff)

    mat = coeff*shc(m)*rad.i4x4lapl(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4lapl(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i4x4lapl2(nr, maxl, m, bc, coeff = 1.0, with_sh_coeff = None):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    mat = coeff*shc(m)*rad.i4x4lapl2(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4lapl2(nr, l, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i2x2coriolis(nr, maxl, m, bc, coeff = 1.0):
    """Create a i2x2 radial operator kronecker with coriolis Q term"""

    cor_r = sh.coriolis_r(maxl, m).tocsr()
    cordr = sh.coriolisdr(maxl, m).tocsr()

    bcr = convert_bc(bc)

    mat = coeff*spsp.kron(cor_r[0,:],rad.i2x1(nr, m, bcr)) + coeff*spsp.kron(cordr[0,:],rad.i2x2d1(nr, m, bcr))
    for ir,l in enumerate(range(m+1, maxl+1)):
        row = coeff*spsp.kron(cor_r[ir+1,:],rad.i2x1(nr, l, bcr)) + coeff*spsp.kron(cordr[ir+1,:],rad.i2x2d1(nr, l, bcr))
        mat = spsp.vstack([mat,row])

    return sphbc.constrain(mat, nr, maxl, m, bc)

def i4x4coriolis(nr, maxl, m, bc, coeff = 1.0):
    """Create a i4x4 radial operator kronecker with coriolis Q term"""

    cor_r = sh.coriolis_r(maxl, m).tocsr()
    cordr = sh.coriolisdr(maxl, m).tocsr()

    bcr = convert_bc(bc)

    mat = coeff*spsp.kron(cor_r[0,:],rad.i4x3(nr, m, bcr)) + coeff*spsp.kron(cordr[0,:],rad.i4x4d1(nr, m, bcr))
    for ir,l in enumerate(range(m+1, maxl+1)):
        row = coeff*spsp.kron(cor_r[ir+1,:],rad.i4x3(nr, l, bcr)) + coeff*spsp.kron(cordr[ir+1,:],rad.i4x4d1(nr, l, bcr))
        mat = spsp.vstack([mat,row])

    return sphbc.constrain(mat, nr, maxl, m, bc)

def qid(nr, maxl, m, qr, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr = convert_bc(bc)

    mat = coeff*rad.qid(nr, m, bcr)
    for l in range(m+1, maxl+1):
        mat = spsp.block_diag((mat,coeff*rad.qid(nr, l, qr, bcr)))

    return sphbc.constrain(mat, nr, maxl, m, bc)
