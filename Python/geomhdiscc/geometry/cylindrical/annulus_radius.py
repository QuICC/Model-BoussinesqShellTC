"""Module provides functions to generate sparse operators for the radial direction in an annulus."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.annulus_radius_boundary as cylbc


def zblk(nr, q, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return cylbc.constrain(mat,bc,q)


def i2x2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 1

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4/(16.0*n*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (a**3*b)/(4.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (a**2*(2.0*b**2*n + a**2 + 2.0*b**2))/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(a**3*b)/(4.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -(a**2*(a**2 + 4.0*b**2))/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n - 1.0)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2 - 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n + 1.0)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 1.0)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return cylbc.constrain(mat, bc, 2)


def i2x2lapl(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(m - n + 2.0)*(m + n - 2.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(2.0*n - 3.0)/(2.0*n)

    # Generate main diagonal
    def d0(n):
        return (a**2*m**2 + a**2*n**2 - 2.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(2.0*n + 3.0)/(2.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(m - n - 2.0)*(m + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return cylbc.constrain(mat, bc, 2)


def i4x4(nr, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-8,9)
    nzrow = 3

    # Generate 8th subdiagonal
    def d_8(n):
        return a**8/(256.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 7th subdiagonal
    def d_7(n):
        return (a**7*b)/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return (3.0*a**6*(2.0*b**2*n + a**2 + 2.0*b**2))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n - 11.0*a**2 - 4.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 - 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 - 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n + 6.0*a**2 + 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**4*(9.0*a**4*n + 33.0*a**4 + 6.0*a**2*b**2*n**2 + 120*a**2*b**2*n + 306.0*a**2*b**2 + 16.0*b**4*n**2 + 80.0*b**4*n + 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a**5*b*(3.0*a**2*n**2 - 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 - 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*(a**4*n**2 - 29.0*a**4 + 16.0*a**2*b**2*n**2 - 304.0*a**2*b**2 + 16.0*b**4*n**2 - 144.0*b**4))/(128.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (a**5*b*(3.0*a**2*n**2 + 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 + 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (a**4*(9.0*a**4*n - 33.0*a**4 - 6.0*a**2*b**2*n**2 + 120.0*a**2*b**2*n - 306.0*a**2*b**2 - 16.0*b**4*n**2 + 80.0*b**4*n - 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n - 6.0*a**2 - 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 + 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 + 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n + 11.0*a**2 + 4.0*b**2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -(3.0*a**6*(a**2 - 2.0*b**2*n + 2.0*b**2))/(64.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 7th superdiagonal
    def d7(n):
        return d_7(n + 3.0)

    # Generate 8th superdiagonal
    def d8(n):
        return d_8(n + 3.0)

    ds = [d_8, d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7, d8]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return cylbc.constrain(mat, bc, 4)


def i4x4lapl(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-6,7)
    nzrow = 3

    # Generate 6th subdiagonal
    def d_6(n):
        return -a**6*(m - n + 6.0)*(m + n - 6.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -a**5*b*(2.0*m**2 - 4.0*n**2 + 41.0*n - 105.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(a**2*m**2*n - 5.0*a**2*m**2 + a**2*n**3 - 3.0*a**2*n**2 - 34.0*a**2*n + 120.0*a**2 - 2.0*b**2*m**2*n - 2.0*b**2*m**2 + 12.0*b**2*n**3 - 90.0*b**2*n**2 + 114.0*b**2*n + 216.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*b*(6.0*a**2*m**2 + 4.0*a**2*n**2 + 25.0*a**2*n - 183.0*a**2 + 16.0*b**2*n**2 - 44.0*b**2*n - 60.0*b**2)/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(a**4*m**2*n + 17.0*a**4*m**2 - a**4*n**3 + 27.0*a**4*n**2 - 8.0*a**4*n - 372.0*a**4 + 16.0*a**2*b**2*m**2*n + 32.0*a**2*b**2*m**2 + 216.0*a**2*b**2*n**2 - 360.0*a**2*b**2*n - 1584.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n - 96.0*b**4)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*b*(2.0*a**2*m**2*n - 16.0*a**2*m**2 + 4.0*a**2*n**3 - 33.0*a**2*n**2 - 55.0*a**2*n + 444.0*a**2 + 8.0*b**2*n**3 - 66.0*b**2*n**2 + 10.0*b**2*n + 348.0*b**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*(a**4*m**2*n**2 - 19.0*a**4*m**2 + a**4*n**4 - 43.0*a**4*n**2 + 396.0*a**4 + 6.0*a**2*b**2*m**2*n**2 - 54.0*a**2*b**2*m**2 + 12.0*a**2*b**2*n**4 - 336.0*a**2*b**2*n**2 + 2052.0*a**2*b**2 + 8.0*b**4*n**4 - 104.0*b**4*n**2 + 288.0*b**4)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a**3*b*(2.0*a**2*m**2*n + 16.0*a**2*m**2 + 4.0*a**2*n**3 + 33.0*a**2*n**2 - 55.0*a**2*n - 444.0*a**2 + 8.0*b**2*n**3 + 66.0*b**2*n**2 + 10.0*b**2*n - 348.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(a**4*m**2*n - 17.0*a**4*m**2 - a**4*n**3 - 27.0*a**4*n**2 - 8.0*a**4*n + 372.0*a**4 + 16.0*a**2*b**2*m**2*n - 32.0*a**2*b**2*m**2 - 216.0*a**2*b**2*n**2 - 360.0*a**2*b**2*n + 1584.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n + 96.0*b**4)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*b*(6.0*a**2*m**2 + 4.0*a**2*n**2 - 25.0*a**2*n - 183.0*a**2 + 16.0*b**2*n**2 + 44.0*b**2*n - 60.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*(a**2*m**2*n + 5.0*a**2*m**2 + a**2*n**3 + 3.0*a**2*n**2 - 34.0*a**2*n - 120.0*a**2 - 2.0*b**2*m**2*n + 2.0*b**2*m**2 + 12.0*b**2*n**3 + 90.0*b**2*n**2 + 114.0*b**2*n - 216.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*b*(2.0*m**2 - 4.0*n**2 - 41.0*n - 105.0)/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -a**6*(m - n - 6.0)*(m + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return cylbc.constrain(mat, bc, 4)


def i4x4lapl2(nr, m, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(m - n + 4.0)*(m - n + 6.0)*(m + n - 6.0)*(m + n - 4.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*b*(2.0*n - 9.0)*(2.0*m**2 - 2.0*n**2 + 18.0*n - 39.0)/(8.0*n*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(a**2*m**4 + 6.0*a**2*m**2*n - 28.0*a**2*m**2 - a**2*n**4 + 10.0*a**2*n**3 - 20.0*a**2*n**2 - 64.0*a**2*n + 192.0*a**2 + 2.0*b**2*m**2*n**2 - 4.0*b**2*m**2*n - 6.0*b**2*m**2 - 6.0*b**2*n**4 + 60.0*b**2*n**3 - 173.0*b**2*n**2 + 46.0*b**2*n + 285.0*b**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(4.0*a**2*m**2*n - 38.0*a**2*m**2 + 12.0*a**2*n**3 - 54.0*a**2*n**2 - 88.0*a**2*n + 461.0*a**2 + 16.0*b**2*n**3 - 72.0*b**2*n**2 + 24.0*b**2*n + 112.0*b**2)/(8.0*n*(n - 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*m**4 + 2.0*a**4*m**2*n**2 - 68.0*a**4*m**2 + 3.0*a**4*n**4 - 68.0*a**4*n**2 + 416.0*a**4 + 8.0*a**2*b**2*m**2*n**2 - 32.0*a**2*b**2*m**2 + 24.0*a**2*b**2*n**4 - 332.0*a**2*b**2*n**2 + 944.0*a**2*b**2 + 8.0*b**4*n**4 - 40.0*b**4*n**2 + 32.0*b**4)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(4.0*a**2*m**2*n + 38.0*a**2*m**2 + 12.0*a**2*n**3 + 54.0*a**2*n**2 - 88.0*a**2*n - 461.0*a**2 + 16.0*b**2*n**3 + 72.0*b**2*n**2 + 24.0*b**2*n - 112.0*b**2)/(8.0*n*(n - 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2*m**4 - 6.0*a**2*m**2*n - 28.0*a**2*m**2 - a**2*n**4 - 10.0*a**2*n**3 - 20.0*a**2*n**2 + 64.0*a**2*n + 192.0*a**2 + 2.0*b**2*m**2*n**2 + 4.0*b**2*m**2*n - 6.0*b**2*m**2 - 6.0*b**2*n**4 - 60.0*b**2*n**3 - 173.0*b**2*n**2 - 46.0*b**2*n + 285.0*b**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*b*(2.0*n + 9.0)*(2.0*m**2 - 2.0*n**2 - 18.0*n - 39.0)/(8.0*n*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*(m - n - 6.0)*(m - n - 4.0)*(m + n + 4.0)*(m + n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return cylbc.constrain(mat, bc, 4)
