"""Module provides functions to generate sparse operators for the radial direction in a spherical shell."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.shell_radius_boundary as radbc


def zblk(nr, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return radbc.constrain(mat,bc)

def x1(nr, a, b, bc, coeff = 1.0, zr = 0):
    """Create operator for x multiplication"""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        return a/2.0

    # Generate diagonal
    def d0(n):
        return b

    # Generate 1st superdiagonal
    def d1(n):
        return a/2.0

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocsr()
    return radbc.constrain(mat, bc)

def x2(nr, a, b, bc, coeff = 1.0, zr = 0):
    """Create operator for x^2 multiplication"""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = -1

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/4.0

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b

    # Generate diagonal
    def d0(n):
        return (a**2 + 2.0*b**2)/2.0

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n)

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocsr()
    return radbc.constrain(mat, bc)

def d1(nr, a, b, bc, coeff = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    row = [2*j for j in range(0,nr)]
    mat = spsp.lil_matrix((nr,nr))
    for i in range(0,nr-1):
        mat[i,i+1:nr:2] = row[i+1:nr:2]
    mat[-zr:,:] = 0

    mat = coeff*(1.0/a)*mat
    return radbc.constrain(mat, bc, location = 'b')

def d2(nr, a, b, bc, coeff = 1.0, zr = 2):
    """Create operator for 2nd derivative"""

    mat = spsp.lil_matrix((nr,nr))
    for i in range(0,nr-2):
        mat[i,i+2:nr:2] = [j*(j**2 - i**2) for j in range(0,nr)][i+2:nr:2]
    mat[-zr:,:] = 0

    mat = coeff*(1.0/a**2)*mat
    return radbc.constrain(mat, bc, location = 'b')

def i1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 1st integral T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return a/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i2(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n*(n - 1.0))

    # Generate diagonal
    def d0(n):
        return -a**2/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n+1.0)

    ds = [d_2, d0, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

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
    return radbc.constrain(mat, bc)

def i2x2lapl(nr, l, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**2*(l - n + 2.0)*(l + n - 1.0))/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a*b*(n - 1.0))/n

    # Generate main diagonal
    def d0(n):
        return (a**2*l**2 + a**2*l + a**2*n**2 - a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(n + 1.0)/n

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(l - n - 1.0)*(l + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

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
    return radbc.constrain(mat, bc)

def i4x4lapl(nr, l, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-6,7)
    nzrow = 3

    # Generate 6th subdiagonal
    def d_6(n):
        return -(a**6*(l - n + 6)*(l + n - 5))/(64*n*(n - 3)*(n - 1)*(n - 2))

    # Generate 5th subdiagonal
    def d_5(n):
        return -(a**5*b*(l**2 + l - 2*n**2 + 19*n - 45))/(16*n*(n - 3)*(n - 1)*(n - 2))

    # Generate 4th subdiagonal
    def d_4(n):
        return (a**4*(a**2*l**2*n - 5*a**2*l**2 + a**2*l*n - 5*a**2*l + a**2*n**3 - 3*a**2*n**2 - 28*a**2*n + 96*a**2 - 2*b**2*l**2*n - 2*b**2*l**2 - 2*b**2*l*n - 2*b**2*l + 12*b**2*n**3 - 84*b**2*n**2 + 96*b**2*n + 192*b**2))/(32*n*(n - 1)*(n - 2)*(n - 3)*(n + 1))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (a**3*b*(3*a**2*l**2 + 3*a**2*l + 2*a**2*n**2 + 11*a**2*n - 75*a**2 + 8*b**2*n**2 - 20*b**2*n - 28*b**2))/(16*n*(n - 1)*(n - 2)*(n + 1))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (a**2*(a**4*l**2*n + 17*a**4*l**2 + a**4*l*n + 17*a**4*l - a**4*n**3 + 24*a**4*n**2 - 5*a**4*n - 294*a**4 + 16*a**2*b**2*l**2*n + 32*a**2*b**2*l**2 + 16*a**2*b**2*l*n + 32*a**2*b**2*l + 192*a**2*b**2*n**2 - 288*a**2*b**2*n - 1344*a**2*b**2 + 16*b**4*n**3 - 112*b**4*n - 96*b**4))/(64*n*(n - 1)*(n - 3)*(n + 2)*(n + 1))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(a**3*b*(a**2*l**2*n - 8*a**2*l**2 + a**2*l*n - 8*a**2*l + 2*a**2*n**3 - 15*a**2*n**2 - 23*a**2*n + 180*a**2 + 4*b**2*n**3 - 30*b**2*n**2 + 2*b**2*n + 156*b**2))/(8*n*(n - 2)*(n - 3)*(n + 2)*(n + 1))

    # Generate main diagonal
    def d0(n):
        return -(a**2*(a**4*l**2*n**2 - 19*a**4*l**2 + a**4*l*n**2 - 19*a**4*l + a**4*n**4 - 37*a**4*n**2 + 312*a**4 + 6*a**2*b**2*l**2*n**2 - 54*a**2*b**2*l**2 + 6*a**2*b**2*l*n**2 - 54*a**2*b**2*l + 12*a**2*b**2*n**4 - 300*a**2*b**2*n**2 + 1728*a**2*b**2 + 8*b**4*n**4 - 104*b**4*n**2 + 288*b**4))/(16*(n - 1)*(n - 2)*(n - 3)*(n + 3)*(n + 2)*(n + 1))

    # Generate 1st superdiagonal
    def d1(n):
        return -(a**3*b*(a**2*l**2*n + 8*a**2*l**2 + a**2*l*n + 8*a**2*l + 2*a**2*n**3 + 15*a**2*n**2 - 23*a**2*n - 180*a**2 + 4*b**2*n**3 + 30*b**2*n**2 + 2*b**2*n - 156*b**2))/(8*n*(n - 1)*(n - 2)*(n + 3)*(n + 2))

    # Generate 2nd superdiagonal
    def d2(n):
        return (a**2*(a**4*l**2*n - 17*a**4*l**2 + a**4*l*n - 17*a**4*l - a**4*n**3 - 24*a**4*n**2 - 5*a**4*n + 294*a**4 + 16*a**2*b**2*l**2*n - 32*a**2*b**2*l**2 + 16*a**2*b**2*l*n - 32*a**2*b**2*l - 192*a**2*b**2*n**2 - 288*a**2*b**2*n + 1344*a**2*b**2 + 16*b**4*n**3 - 112*b**4*n + 96*b**4))/(64*n*(n - 1)*(n - 2)*(n + 3)*(n + 1))

    # Generate 3rd superdiagonal
    def d3(n):
        return (a**3*b*(3*a**2*l**2 + 3*a**2*l + 2*a**2*n**2 - 11*a**2*n - 75*a**2 + 8*b**2*n**2 + 20*b**2*n - 28*b**2))/(16*n*(n - 1)*(n + 2)*(n + 1))

    # Generate 4th superdiagonal
    def d4(n):
        return (a**4*(a**2*l**2*n + 5*a**2*l**2 + a**2*l*n + 5*a**2*l + a**2*n**3 + 3*a**2*n**2 - 28*a**2*n - 96*a**2 - 2*b**2*l**2*n + 2*b**2*l**2 - 2*b**2*l*n + 2*b**2*l + 12*b**2*n**3 + 84*b**2*n**2 + 96*b**2*n - 192*b**2))/(32*n*(n - 1)*(n + 3)*(n + 2)*(n + 1))

    # Generate 5th superdiagonal
    def d5(n):
        return -(a**5*b*(l**2 + l - 2*n**2 - 19*n - 45))/(16*n*(n + 3)*(n + 2)*(n + 1))

    # Generate 6th superdiagonal
    def d6(n):
        return -(a**6*(l - n - 5)*(l + n + 6))/(64*n*(n + 3)*(n + 2)*(n + 1))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i4x4lapl2(nr, l, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian^2 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return (a**4*(l - n + 6.0)*(l + n - 5.0)*(l - n + 4.0)*(l + n - 3.0))/(16.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(a**3*b*(n - 4.0)*(l**2 + l - n**2 + 8.0*n - 15.0))/(2.0*n*(n - 1.0)*(n - 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**2*(a**2*l**4 + 2.0*a**2*l**3 + 5.0*a**2*l**2*n - 20.0*a**2*l**2 + 5.0*a**2*l*n - 21.0*a**2*l - a**2*n**4 + 9.0*a**2*n**3 - 17.0*a**2*n**2 - 39.0*a**2*n + 108.0*a**2 + 2.0*b**2*l**2*n**2 - 4.0*b**2*l**2*n - 6.0*b**2*l**2 + 2.0*b**2*l*n**2 - 4.0*b**2*l*n - 6.0*b**2*l - 6.0*b**2*n**4 + 54.0*b**2*n**3 - 138.0*b**2*n**2 + 18.0*b**2*n + 216.0*b**2))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a*b*(a**2*l**2*n - 8*a**2*l**2 + a**2*l*n - 8*a**2*l + 3*a**2*n**3 - 12*a**2*n**2 - 15*a**2*n + 72*a**2 + 4*b**2*n**3 - 16*b**2*n**2 + 4*b**2*n + 24*b**2))/(2*n*(n + 1)*(n - 2))

    # Generate main diagonal
    def d0(n):
        return (3*a**4*l**4 + 6*a**4*l**3 + 2*a**4*l**2*n**2 - 47*a**4*l**2 + 2*a**4*l*n**2 - 50*a**4*l + 3*a**4*n**4 - 51*a**4*n**2 + 228*a**4 + 8*a**2*b**2*l**2*n**2 - 32*a**2*b**2*l**2 + 8*a**2*b**2*l*n**2 - 32*a**2*b**2*l + 24*a**2*b**2*n**4 - 264*a**2*b**2*n**2 + 672*a**2*b**2 + 8*b**4*n**4 - 40*b**4*n**2 + 32*b**4)/(8*(n - 1)*(n - 2)*(n + 2)*(n + 1))

    # Generate 1st superdiagonal
    def d1(n):
        return (a*b*(a**2*l**2*n + 8*a**2*l**2 + a**2*l*n + 8*a**2*l + 3*a**2*n**3 + 12*a**2*n**2 - 15*a**2*n - 72*a**2 + 4*b**2*n**3 + 16*b**2*n**2 + 4*b**2*n - 24*b**2))/(2*n*(n + 2)*(n - 1))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(a**2*(a**2*l**4 + 2*a**2*l**3 - 5*a**2*l**2*n - 20*a**2*l**2 - 5*a**2*l*n - 21*a**2*l - a**2*n**4 - 9*a**2*n**3 - 17*a**2*n**2 + 39*a**2*n + 108*a**2 + 2*b**2*l**2*n**2 + 4*b**2*l**2*n - 6*b**2*l**2 + 2*b**2*l*n**2 + 4*b**2*l*n - 6*b**2*l - 6*b**2*n**4 - 54*b**2*n**3 - 138*b**2*n**2 - 18*b**2*n + 216*b**2))/(4*n*(n + 3)*(n - 1)*(n + 1))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(a**3*b*(n + 4)*(l**2 + l - n**2 - 8*n - 15))/(2*n*(n + 2)*(n + 1))

    # Generate 4th superdiagonal
    def d4(n):
        return (a**4*(l + n + 6)*(l - n - 3)*(l + n + 4)*(l - n - 5))/(16*n*(n + 3)*(n + 2)*(n + 1))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i2x1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of x T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -d_3(n + 1.0)

    # Generate main diagonal
    def d0(n):
        return -a**2*b/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_3(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n + 1.0)

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i2x2d1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 D_x T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(n - 2.0)/(2.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*n + 3.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return a**2*b/((n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**2*n - 3.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(n + 2.0)/(2.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(n + 3.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i4x3(nr, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^3 T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-7,8)
    nzrow = 3

    # Generate 7th subdiagonal
    def d_7(n):
        return a**7/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return 3.0*a**6*b/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -a**5*(a**2*n - 11.0*a**2 - 12.0*b**2*n - 12.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -a**4*b*(3.0*a**2*n - 15.0*a**2 - 2.0*b**2*n - 2.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -3.0*a**5*(a**2*n + 6.0*a**2 + 12.0*b**2*n + 24.0*b**2)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(3.0*a**2*n + 51.0*a**2 + 16.0*b**2*n + 32.0*b**2)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 3.0*a**5*(a**2*n**2 - 5.0*a**2*n - 34.0*a**2 + 8.0*b**2*n**2 - 40.0*b**2*n - 192.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return 3.0*a**4*b*(a**2*n**2 - 19.0*a**2 + 2.0*b**2*n**2 - 18.0*b**2)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**5*(a**2*n**2 + 5.0*a**2*n - 34.0*a**2 + 8.0*b**2*n**2 + 40.0*b**2*n - 192.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4*b*(3.0*a**2*n - 51.0*a**2 + 16.0*b**2*n - 32.0*b**2)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3.0*a**5*(a**2*n - 6.0*a**2 + 12.0*b**2*n - 24.0*b**2)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*b*(3.0*a**2*n + 15.0*a**2 - 2.0*b**2*n + 2.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(a**2*n + 11.0*a**2 - 12.0*b**2*n + 12.0*b**2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return 3.0*a**6*b/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d7(n):
        return a**7/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def i4x4d1(nr, a, b, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 D_x T_n(x)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-7,8)
    nzrow = 3

    # Generate 7th subdiagonal
    def d_7(n):
        return a**7*(n - 7.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return a**6*b*(n - 6.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(n - 5.0)*(a**2*n + 13.0*a**2 + 24.0*b**2*n + 24.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b*(n - 4.0)*(3.0*a**2 + b**2*n + b**2)/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*(3.0*a**4*n**2 - 15.0*a**4*n - 102.0*a**4 + 24.0*a**2*b**2*n**2 - 216.0*a**2*b**2*n - 528.0*a**2*b**2 - 16.0*b**4*n**2 - 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(3.0*a**2*n**2 - 3.0*a**2*n - 78.0*a**2 + 8.0*b**2*n**2 - 24.0*b**2*n - 80.0*b**2)/(16.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*a**3*(a**4*n**3 + 10.0*a**4*n**2 - 29.0*a**4*n - 190.0*a**4 + 16.0*a**2*b**2*n**3 + 96.0*a**2*b**2*n**2 - 304.0*a**2*b**2*n - 1344.0*a**2*b**2 + 16.0*b**4*n**3 + 32.0*b**4*n**2 - 144.0*b**4*n - 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return -3.0*a**4*b*(a**2*n**2 - 14.0*a**2 + 2.0*b**2*n**2 - 18.0*b**2)/(2*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**3*(a**4*n**3 - 10.0*a**4*n**2 - 29.0*a**4*n + 190.0*a**4 + 16.0*a**2*b**2*n**3 - 96.0*a**2*b**2*n**2 - 304.0*a**2*b**2*n + 1344.0*a**2*b**2 + 16.0*b**4*n**3 - 32.0*b**4*n**2 - 144.0*b**4*n + 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**4*b*(3.0*a**2*n**2 + 3.0*a**2*n - 78.0*a**2 + 8.0*b**2*n**2 + 24.0*b**2*n - 80.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(3.0*a**4*n**2 + 15.0*a**4*n - 102.0*a**4 + 24.0*a**2*b**2*n**2 + 216.0*a**2*b**2*n - 528.0*a**2*b**2 - 16.0*b**4*n**2 + 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*b*(n + 4.0)*(3.0*a**2 - b**2*n + b**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(n + 5.0)*(a**2*n - 13.0*a**2 + 24.0*b**2*n - 24.0*b**2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -a**6*b*(n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d7(n):
        return -a**7*(n + 7.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def qid(nr, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    offsets = [0]
    diags = [[0]*q + [1]*(nr-q)]

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, bc)

def linear_r2x(ro, rratio):
    """Calculat a and b for linear map r = a*x + b"""

    b = (ro*rratio + ro)/2.0;
    a = ro - b;

    return (a, b)
