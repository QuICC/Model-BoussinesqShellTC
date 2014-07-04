"""Module provides functions to generate sparse operators for the radial direction in a sphere."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.sphere_radius_boundary as sphbc


def zblk(nr, l, q, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return sphbc.constrain(mat,l,bc,q)


def i2x2(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(16.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -1.0/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 1)


def i2x2lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 Laplacian T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return  -(l - n + 2.0)*(l + n - 1.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (l**2 + l + n**2 - 1.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(l - n - 1.0)*(l + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 1)


def i4x4(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return 1.0/(256.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return 3.0/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(n**2 - 19.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(3.0*(3.0*n + 11.0))/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*(n**2 - 29.0))/(128.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (3.0*(3.0*n - 11.0))/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n**2 - 19.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -d_3(n + 2.0)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 3.0)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 2)


def i4x4lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(l - n + 6.0)*(l + n - 5.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (l**2*n - 5.0*l**2 + l*n - 5.0*l + n**3 - 3.0*n**2 - 28.0*n + 96.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (l**2*n + 17.0*l**2 + l*n + 17.0*l - n**3 + 24.0*n**2 - 5.0*n - 294.0)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -(l**2*n**2 - 19.0*l**2 + l*n**2 - 19.0*l + n**4 - 37.0*n**2 + 312.0)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (l**2*n - 17.0*l**2 + l*n - 17.0*l - n**3 - 24.0*n**2 - 5.0*n + 294.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (l**2*n + 5.0*l**2 + l*n + 5.0*l + n**3 + 3.0*n**2 - 28.0*n - 96.0)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(l - n - 5.0)*(l + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 2)


def i4x4lapl2(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian^2 T_n(x)."""

    parity = l%2
    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return ((l - n + 4.0)*(l - n + 6.0)*(l + n - 5.0)*(l + n - 3.0))/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -((l - n + 4.0)*(l + n - 3.0)*(l**2 + l + n**2 - 2.0*n - 9.0))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*l**4 + 6.0*l**3 + 2.0*l**2*n**2 - 47.0*l**2 + 2.0*l*n**2 - 50.0*l + 3.0*n**4 - 51.0*n**2 + 228.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -((l - n - 3.0)*(l + n + 4.0)*(l**2 + l + n**2 + 2.0*n - 9.0))/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return ((l - n - 5.0)*(l - n - 3.0)*(l + n + 4.0)*(l + n + 6.0))/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 2)


def i2x1(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of x T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -d_2(n + 1.0)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_2(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 1)


def i2x2d1(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 D T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return (n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (n + 3.0)/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_2(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -d_1(n)

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 1)


def i4x3(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of x^3 T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-3,5)
    else:
        offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return 1.0/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(n - 11.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -3.0*(n + 6.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 3.0*(n**2 - 5.0*n - 34.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*(n**2 + 5.0*n - 34.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -3.0*(n - 6.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(n + 11.0)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return 1.0/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 2)


def i4x4d1(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 D T_n(x)."""

    ns = np.arange((l+1)%2, 2*nr, 2)
    if l%2 == 0:
        offsets = np.arange(-3,5)
    else:
        offsets = np.arange(-4,4)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return (n - 7.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (n - 5.0)*(n + 13.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -3.0*(n**2 - 5.00*n - 34.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*(n**3 + 10.0*n**2 - 29.0*n - 190.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*(n**3 - 10.0*n**2 - 29.0*n + 190.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 3.0*(n**2 + 5.0*n - 34.0)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(n - 13.0)*(n + 5.0)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(n + 7.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (l+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, 2)


def qid(nr, l, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    offsets = [0]
    diags = [[0]*q + [1]*(nr-q)]

    mat = coeff*spsp.diags(diags, offsets)
    return sphbc.constrain(mat, l, bc, q)