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
    mat = coeff*spsp.kron(c1d.i1d1(nz,bcz), rad.i1x1(nr, a, b, bcr))
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
    """Create a i2x2d1 radial operator kronecker with an i2"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x2d1(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2e1(nr, nz, a, b, bc, coeff = 1.0, zscale = 1.0):
    """Create a i2x2 in R kronecker with an i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(nz,bcz), rad.i2x2(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2lapl(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    print("NOT RIGHT")
    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), rad.i2x2lapl(nr, m, a, b, bcr))
    mat = mat + coeff*spsp.kron(c1d.i2d2(nz,bcz), rad.i2x2(nr, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i4j4x4(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i4x4 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4x4(nr, m, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j4x4lapl(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    print("NOT RIGHT")
    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4x4lapl(nr, m, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j4x4lapl2(nr, nz, m, a, b, bc, coeff = 1.0):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    print("NOT RIGHT")
    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), rad.i4x4lapl2(nr, m, a, b, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def qid(nr, nz, m, qr, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), rad.qid(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, qr, qz, bc)
