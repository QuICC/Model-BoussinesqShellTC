"""Module provides functions to generate sparse operators in a spherical shell using a spherical harmonics expansion in the angular directions."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp
import geomhdiscc.geometry.spherical.shell_radius as rad
import geomhdiscc.geometry.spherical.shell_boundary as sphbc


def zblk(nr, nl, qr, bc):
    """Create a block of zeros"""

    mat = spsp.kron(rad.zblk(nl,0,[0]),rad.zblk(nr,qr,bc))
    return sphbc.constrain(mat, nr, nl, bc)


def i2x2(nr, nl, a, b, bc, coeff = 1.0):
    """Create a i2x2 radial operator kronecker with an identity"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.i2x2(nr, a, b, bc))
    return sphbc.constrain(mat, nr, nl, bc)


def i2x2lapl(nr, nl, l, a, b, bc, coeff = 1.0):
    """Create a i2x2lapl radial operator kronecker with an identity"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.i2x2lapl(nr, l, a, b, bc))
    return sphbc.constrain(mat, nr, nl, bc)


def i4x4(nr, nl, a, b, bc, coeff = 1.0):
    """Create a i4x4 radial operator kronecker with an identity"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.i4x4(nr, a, b, bc))
    return sphbc.constrain(mat, nr, nl, bc)


def i4x4lapl(nr, nl, l, a, b, bc, coeff = 1.0):
    """Create a i4x4lapl radial operator kronecker with an identity"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.i4x4lapl(nr, l, a, b, bc))
    return sphbc.constrain(mat, nr, nl, bc)


def i4x4lapl2(nr, nl, l, a, b, bc, coeff = 1.0):
    """Create a i4x4lapl2 radial operator kronecker with an identity"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.i4x4lapl2(nr, l, a, b, bc))
    return sphbc.constrain(mat, nr, nl, bc)


def qid(nr, nl, qr, bc, coeff = 1.0):
    """Create a quasi identity block order qr in r"""

    mat = coeff*spsp.kron(rad.qid(nl,0,[0]), rad.qid(nr,qr,bc))
    return sphbc.constrain(mat, nr, nl, bc)
