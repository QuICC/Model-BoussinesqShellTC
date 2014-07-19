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

def zblk(nr, nz, bc):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,qz,bcz),rad.zblk(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, bc)

def i2j2x2(nr, nz, bc, coeff = 1.0):
    """Create a i2x2 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.j2(nz,2,bcz), rad.i2x2(nr, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i2j2x2lapl(nr, nz, m, bc, coeff = 1.0):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz), rad.i2x2lapl(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, 2, 2, bc)

def i4j4x4(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,4,bcz), rad.i4x4(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j4x4lapl(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,4,bcz), rad.i4x4lapl(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def i4j4x4lapl2(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,4,bcz), rad.i4x4lapl2(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, 4, 4, bc)

def qid(nr, nz, m, qr, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), rad.qid(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, qr, qz, bc)
