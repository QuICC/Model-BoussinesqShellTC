"""Module provides functions to generate sparse operators in a cylinder."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cylindrical.cylinder_radius as rad
import geomhdiscc.geometry.cylindrical.cylinder_boundary as cylbc


def convert_bc(bc):
    """Convert boundary dictionary into x and z kronecker product boundaries"""

    if bc['r']bcz < 0:
        bcx = bc['r']
    else:
        bcx = bcz

    if bc['z']bcz < 0:
        bcz = bc['z']
    else:
        bcz = bcz

    return (bcr, bcz)


def zblk(nr, nz, qr, qz, bc):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,qz,bcz),rad.zblk(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, bc, qr, qz)


def i2j2x2(nr, nz, m, bc, coeff = 1.0):
    """Create a i2x2 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.j2(nz,2,bcz), rad.i2x2(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, bc, 1, 2)


def i2j2x2lapl(nr, nz, m, bc, coeff = 1.0):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz), rad.i2x2lapl(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, bc, 1, 2)


def i4j4x4(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,4,bcz), rad.i4x4(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, bc, 2, 4)


def i4j4x4lapl(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,4,bcz), rad.i4x4lapl(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, bc, 2, 4)


def i4j4x4lapl2(nr, nz, m, bc, coeff = 1.0):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,4,bcz), rad.i4x4lapl2(nr, m, bcr))
    return cylbc.constrain(mat, nr, nz, bc, 2, 4)


def qid(nr, nz, m, qr, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,0,bcz), rad.qid(nr,qr,bcr))
    return cylbc.constrain(mat, nr, nz, bc, qr, qz)
