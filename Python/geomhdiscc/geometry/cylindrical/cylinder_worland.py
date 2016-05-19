"""Module provides functions to generate sparse operators in a cylinder with Worland expansion in radius and Chebyshev in the vertical."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cylindrical.cylinder_radius_worland as rad
import geomhdiscc.geometry.cylindrical.cylinder_boundary_worland as cylbc


def convert_bc(bc):
    """Convert boundary dictionary into r and z kronecker product boundaries"""

    if bc['r'][0] < 0:
        bcr = bc['r']
    else:
        bcr = rad.radbc.no_bc()

        for key, val in bc['r'].items():
            if key != 0:
                bcr[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcr, bcz)

def zblk(nr, nz, m, qr, qz, bc, restriction = None):
    """Create a block of zeros"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.zblk(nz,bcz),rad.zblk(nr,m,bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, qr, qz, bc, restriction = restriction)

def i2j2(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i2(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 1, 2, bc, restriction = restriction)

def i2laplhj2(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i2laph in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i2laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 1, 2, bc, restriction = restriction)

def i2j2lapl(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i2 in R kronecker with i2 in Z of the laplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i2laplh(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i2d2(nz, bcz, cscale = zscale), rad.i2(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, nz, m, 1, 2, bc, restriction = restriction)

def i4j2(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i4 in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i4(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 2, 2, bc, restriction = restriction)

def i4laplhj2(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i4laplh in R kronecker with i2 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i4laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 2, 2, bc, restriction = restriction)

def i4laplhj2e1(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i4laplh in R kronecker with i2d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i2d1(nz, bcz, cscale = zscale), rad.i4laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 2, 2, bc, restriction = restriction)

def i4j2lapl2(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i4 in R kronecker with i2 in Z of the bilaplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i2(nz, bcz), rad.i4lapl2h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i2d2(nz, bcz, cscale = zscale), rad.i4laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, nz, m, 2, 2, bc, restriction = restriction)

def i6j4(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, bcz), rad.i6(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 3, 4, bc, restriction = restriction)

def i6laplhj4(nr, nz, m, bc, coeff = 1.0, restriction = None):
    """Create a i6laplh in R kronecker with i4 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4(nz, bcz), rad.i6laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 3, 4, bc, restriction = restriction)

def i6laplhj4e1(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i6laplh in R kronecker with i4d1 in Z"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.i4d1(nz, bcz, cscale = zscale), rad.i6laplh(nr, m, bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, 3, 4, bc, restriction = restriction)

def i6j4lapl2(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z of the bilaplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, bcz), rad.i6lapl2h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, bcz, cscale = zscale), rad.i6laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, nz, m, 3, 4, bc, restriction = restriction)

def i6j4lapl3(nr, nz, m, bc, coeff = 1.0, zscale = 1.0, restriction = None):
    """Create a i6 in R kronecker with i4 in Z of the trilaplacian"""

    bcr, bcz = convert_bc(bc)
    mat = utils.restricted_kron_2d(c1d.i4(nz, bcz), rad.i6lapl3h(nr, m, bcr), restriction = restriction)
    bcr[0] = min(bcr[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + utils.restricted_kron_2d(c1d.i4d2(nz, bcz, 2.0, cscale = zscale), rad.i6lapl2h(nr, m, bcr), restriction = restriction)
    mat = mat + utils.restricted_kron_2d(c1d.i4d4(nz, bcz, cscale = zscale), rad.i6laplh(nr, m, bcr), restriction = restriction)
    mat *= coeff
    return cylbc.constrain(mat, nr, nz, m, 3, 4, bc, restriction = restriction)

def qid(nr, nz, m, qr, qz, bc, coeff = 1.0, restriction = None):
    """Create a quasi identity block order qr in r"""

    bcr, bcz = convert_bc(bc)
    mat = coeff*utils.restricted_kron_2d(c1d.qid(nz,qz,bcz), rad.qid(nr, m, qr,bcr), restriction = restriction)
    return cylbc.constrain(mat, nr, nz, m, qr, qz, bc, restriction = restriction)
