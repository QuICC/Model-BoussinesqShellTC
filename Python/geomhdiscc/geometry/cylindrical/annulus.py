"""Module provides functions to generate sparse operators in a cylindrical annulus."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cylindrical.annulus_radius as rad
import geomhdiscc.geometry.cylindrical.annulus_boundary as cylbc


def convert_bc(bc):
    """Convert boundary dictionary into r and z kronecker product boundaries"""

    if bc['r'][0] < 0:
        bcr = bc['r']
    else:
        bcr = rad.radbc.no_bc()
        bcr['r'] = bc['r'].get('r',0)

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        bcz['r'] = bc['z'].get('r',0)

    return (bcr, bcz)

def zblk(nr, nz, qr, qz, bc):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,bcz),rad.zblk(nr,bcr))
    return cylbc.constrain(mat, nr, nz, qr, qz, bc)

def i1j1(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i1 in R kronecker with an i1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), rad.i1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i1j1x1d1(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i1x1d1 in R kronecker with an i1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), rad.i1x1d1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i1j1x1e1(nr, nz, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i1x1 in R kronecker with an i1d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = zscale*coeff*spsp.kron(c1d.i1d1(nz,bcz), rad.i1x1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i2 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x1(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i2x1 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i2x2 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x2(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2d1(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i2x2d1 in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x2d1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2e1(nr, nz, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2x2 in R kronecker with an i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*zscale*spsp.kron(c1d.i2d1(nz,bcz), rad.i2x2(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2div(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i2x2div in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x2div(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j0x2laplh(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i2x2laplh in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz, 0, bcz), rad.i2x2laplh(nr, m, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 0, bc)

def i2j2x2lapl(nr, nz, m, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2x2lapl in R kronecker with an i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), rad.i2x2laplh(nr, m, a, b, bcr))
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + zscale**2*spsp.kron(c1d.i2d2(nz,bcz), rad.i2x2(nr, a, b, bcr))
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i4j4(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i4 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j4x4(nr, nz, a, b, bc, coeff = 1.0):
    """Create a i4x4 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4x4(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j0x4laplh(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i4x4laplh in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz, 0, bcz), rad.i4x4laplh(nr, m, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 0, bc)

def i4j4x4lapl(nr, nz, m, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i4x4lapl in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), rad.i4x4laplh(nr, m, a, b, bcr))
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + zscale**2*spsp.kron(c1d.i4d2(nz,bcz), rad.i4x4(nr, a, b, bcr))
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j0x4lapl2h(nr, nz, m, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i4x4lapl2h in R kronecker with an identity in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.qid(nz, 0, bcz), rad.i4x4lapl2h(nr, m, a, b, bcr))
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, 4, 0, bc)

def i4j4x4lapl2(nr, nz, m, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i4x4lapl2 in R kronecker with an i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), rad.i4x4lapl2h(nr, m, a, b, bcr))
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + 2.0*zscale**2*spsp.kron(c1d.i4d2(nz,bcz), rad.i4x4laplh(nr, m, a, b, bcr))
    mat = mat + zscale**4*spsp.kron(c1d.i4d4(nz,bcz), rad.i4x4(nr, a, b, bcr))
    mat = coeff*mat
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def qid(nr, nz, m, qr, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), rad.qid(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, qr, qz, bc)
