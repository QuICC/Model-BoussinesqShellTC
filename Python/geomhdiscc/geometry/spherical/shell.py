"""Module provides functions to generate sparse operators in a spherical shell using a spherical harmonics expansion in the angular directions."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.geometry.spherical.shell_radius as rad
import geomhdiscc.geometry.spherical.shell_sh as sh
import geomhdiscc.geometry.spherical.shell_boundary as sphbc


def convert_bc(bc):
    """Convert boundary condition to be suitable for kronecker product boundaries"""

    if bc[0] < 0:
        bcr = bc
    else:
        bcr = rad.radbc.no_bc()
        for k, v in bc.items():
            if k != 0:
                bcr[k] = v

    return bcr

def sh_coeff(coeff):
    """Compute coefficients from spherical harmonic expansion"""

    if coeff == 'laplh':
        def fct(x):
            return x*(x + 1.0)
    elif coeff == 'laplh_1':
        def fct(x):
            return 1.0/(x*(x + 1.0))
    else:
        def fct(x):
            return 1.0
    
    return fct

def fix_l_zero(nr, m, mat, bc, fix):
    """Fix problems with unused l = 0 modes"""

    if m > 0 or not fix:
        return mat
    elif fix == 'zero':
        return rad.zblk(nr, bc)
    elif fix == 'set':
        return rad.qid(nr, 0, bc)
    else:
        raise RuntimeError("Unkown l=0 fix!")

def make_sh_operator(op, nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False)
    """Generic function to create a coupled spherical harmonics operator"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*op(nr, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*op(nr, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def make_sh_loperator(op, nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False)
    """Generic function to create a coupled l dependent spherical harmonics operator"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*op(nr, m, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*op(nr, l, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def zblk(nr, maxnl, m, bc):
    """Create a block of zeros"""

    bcr = convert_bc(bc)

    nl = maxnl - m
    mat = spsp.kron(rad.zblk(nl,rad.radbc.no_bc()),rad.zblk(nr,bcr))
    return sphbc.constrain(mat, nr, maxnl, m, bc)

def i2(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i2(nr, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2(nr, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i2x2(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2x2 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i2x2(nr, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x2(nr, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i2x3(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2x3 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i2x3(nr, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x3(nr, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i2x2coriolis(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2x2 radial operator kronecker with coriolis Q term"""

    # Only compute if there are at least 2 harmonic degrees
    if maxnl - m > 1:
        cor_r = sh.coriolis_r(maxnl, m).tocsr()
        cordr = sh.coriolisdr(maxnl, m).tocsr()

        bcr = convert_bc(bc)
        shc = sh_coeff(with_sh_coeff)

        assert(l_zero_fix == 'zero')
        bcr = sphbc.ldependent_bc(bcr, m)
        rmat1 = rad.i2x1(nr, a, b, bcr)
        rmat2 = rad.i2x2d1(nr, a, b, bcr)
        rmat1 = fix_l_zero(nr, m, rmat1, bcr, l_zero_fix)
        rmat2 = fix_l_zero(nr, m, rmat2, bcr, l_zero_fix)
        mat = coeff*shc(m)*spsp.kron(cor_r[0,:],rmat1) + coeff*shc(m)*spsp.kron(cordr[0,:], rmat2)
        for ir,l in enumerate(range(m+1, maxnl)):
            bcr = sphbc.ldependent_bc(bcr, l)
            row = coeff*shc(l)*spsp.kron(cor_r[ir+1,:],rad.i2x1(nr, a, b, bcr)) + coeff*shc(l)*spsp.kron(cordr[ir+1,:],rad.i2x2d1(nr, a, b, bcr))
            mat = spsp.vstack([mat,row])
    else:
        mat = rad.zblk(nr, bc)

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i2x2lapl(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i2x2lapl(nr, m, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x2lapl(nr, l, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i2x3lapl(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i2x3lapl radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i2x3lapl(nr, m, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i2x3lapl(nr, l, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i4x4(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i4x4 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i4x4(nr, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4(nr, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i4x4coriolis(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i4x4 radial operator kronecker with coriolis Q term"""

    # Only compute if there are at least 2 harmonic degrees
    if maxnl - m > 1:
        cor_r = sh.coriolis_r(maxnl, m).tocsr()
        cordr = sh.coriolisdr(maxnl, m).tocsr()

        bcr = convert_bc(bc)
        shc = sh_coeff(with_sh_coeff)

        assert(l_zero_fix == 'zero')
        bcr = sphbc.ldependent_bc(bcr, m)
        rmat1 = rad.i4x3(nr, a, b, bcr)
        rmat2 = rad.i4x4d1(nr, a, b, bcr)
        rmat1 = fix_l_zero(nr, m, rmat1, bcr, l_zero_fix)
        rmat2 = fix_l_zero(nr, m, rmat2, bcr, l_zero_fix)
        mat = coeff*shc(m)*spsp.kron(cor_r[0,:],rmat1) + coeff*shc(m)*spsp.kron(cordr[0,:],rmat2)
        for ir,l in enumerate(range(m+1, maxnl)):
            bcr = sphbc.ldependent_bc(bcr, l)
            row = coeff*shc(l)*spsp.kron(cor_r[ir+1,:],rad.i4x3(nr, a, b, bcr)) + coeff*shc(l)*spsp.kron(cordr[ir+1,:],rad.i4x4d1(nr, a, b, bcr))
            mat = spsp.vstack([mat,row])

    else:
        mat = rad.zblk(nr, bc)

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i4x4lapl(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i4x4lapl(nr, m, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4lapl(nr, l, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def i4x4lapl2(nr, maxnl, m, a, b, bc, coeff = 1.0, with_sh_coeff = None, l_zero_fix = False):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    bcr = convert_bc(bc)
    shc = sh_coeff(with_sh_coeff)

    bcr = sphbc.ldependent_bc(bcr, m)
    mat = coeff*shc(m)*rad.i4x4lapl2(nr, m, a, b, bcr)
    mat = fix_l_zero(nr, m, mat, bcr, l_zero_fix)
    for l in range(m+1, maxnl):
        bcr = sphbc.ldependent_bc(bcr, l)
        mat = spsp.block_diag((mat,coeff*shc(l)*rad.i4x4lapl2(nr, l, a, b, bcr)))

    return sphbc.constrain(mat, nr, maxnl, m, bc, l_zero_fix)

def qid(nr, maxnl, m, qr, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr = convert_bc(bc)

    nl = maxnl - m
    mat = coeff*spsp.kron(rad.qid(nl,0,rad.radbc.no_bc()), rad.qid(nr,qr,bcr))
    return sphbc.constrain(mat, nr, maxnl, m, bc)
