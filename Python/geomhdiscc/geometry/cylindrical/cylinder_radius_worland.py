"""Module provides functions to generate sparse operators for the radial direction in a cylinder with Worland expansion."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.cylinder_radius_boundary_worland as radbc
import geomhdiscc.geometry.worland.worland_basis as wb


def zblk(nr, m, bc):
    """Create a block of zeros"""

    # Copy BC dict as we modify it!
    bc = dict(bc)

    mat = spsp.coo_matrix((nr,nr))
    return radbc.constrain(mat,m,bc)

def i2(nr, m, bc, coeff = 1.0):
    """Create operator for 2nd integral r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    # Copy BC dict as we modify it!
    bc = dict(bc)

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        if m == 0:
            return wb.worland_norm_row(n,0,-2)/((2.0*n - 3.0)*(2.0*n - 1.0))
        else:
            return wb.worland_norm_row(n,m,-2)*4.0*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -wb.worland_norm_row(n,m,-1)*8.0*m*(m + n - 1.0)/((m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n,m,0)*2.0*(2.0*m**2 - 4.0*m*n - 4.0*n**2 + 1.0)/((m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n,m,1)*2.0*m*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)/((m + n)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n,m,2)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)/(4.0*(m + n)*(m + n + 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, m, bc)

def i2r_1dr(nr, m, bc, coeff = 1.0):
    """Create operator for 2nd integral of 1/r D r r P_n^{-1/2,1/2}(2r^2-1) projected onto P_n^{-1/2,-1/2}(2r^2-1)."""
    if m != 0:
        raise RuntimeError("Operator only exists for m = 0!")

    # Copy BC dict as we modify it!
    bc = dict(bc)

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n,m,-2)*4.0*(n - 1.0)/((2.0*n - 3.0)*(2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n,m,-1)*2.0/(2*n - 1.0)

    # Generate diagonal
    def d0(n):
        return -wb.worland_norm_row(n,m,0)*1.0/n

    # Generate 1st superdiagonal
    def d1(n):
        return -wb.worland_norm_row(n,m,1)*(2.0*n + 1.0)/(2.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, m+1, bc)

def i2dr(nr, m, bc, coeff = 1.0):
    """Create operator for 2nd integral of D P_n^{-1/2,-1/2}(2r^2-1) projected onto r P_n^{-1/2,1/2}(2r^2-1)."""
    if m != 1:
        raise RuntimeError("Operator only exists for m = 1!")

    # Copy BC dict as we modify it!
    bc = dict(bc)

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,3)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n,m,-1)*2.0/(2.0*n - 1.0)

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n,m,0)*1.0/(n + 1.0)

    # Generate 1st superdiagonal
    def d1(n):
        return -wb.worland_norm_row(n,m,1)*(2.0*n + 1.0)/(2.0*n*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -wb.worland_norm_row(n,m,2)*(2.0*n + 1.0)*(2.0*n + 3.0)/(4.0*(n + 1.0)**2*(n + 2.0))

    ds = [d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, m-1, bc)

def i2laplh(nr, m, bc, coeff = 1.0):
    """Create operator for 2nd integral of laplh r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n,m,-1)*16.0*(m + n - 1.0)**2/((m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n,m,0)*8.0*(4.0*m*n + m + 4.0*n**2 - 1.0)/((m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n,m,1)*4.0*(n + 1.0)*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)/((m + n)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, m, bc)

def i4(nr, m, bc, coeff = 1.0):
    """Create operator for 4th integral r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        if m == 0:
            return wb.worland_norm_row(n, 0, -4)/((2.0*n - 7.0)*(2.0*n - 5.0)*(2.0*n - 3.0)*(2.0*n - 1.0))
        else:
            return wb.worland_norm_row(n, m, -4)*16.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -wb.worland_norm_row(n, m, -3)*64.0*m*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*16.0*(m + n - 2.0)*(m + n - 1.0)*(6.0*m**2 - 4.0*m*n + 4.0*m - 4.0*n**2 + 8.0*n + 5.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -wb.worland_norm_row(n, m, -1)*16.0*m*(m + n - 1.0)*(4.0*m**2 - 12.0*m*n + 6.0*m - 12.0*n**2 + 12.0*n + 17.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*2.0*(8.0*m**4 - 96.0*m**3*n - 48.0*m**2*n**2 + 100.0*m**2 + 96.0*m*n**3 - 120.0*m*n + 48.0*n**4 - 120.0*n**2 + 27.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*4.0*m*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(4.0*m**2 - 12.0*m*n - 6.0*m - 12.0*n**2 - 12.0*n + 17.0)/((m + n)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(6.0*m**2 - 4.0*m*n - 4.0*m - 4.0*n**2 - 8.0*n + 5.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 4th superdiagonal
    def d4(n):
        return wb.worland_norm_row(n, m, 4)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)/(16.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, m, bc)

def i4dr(nr, m, bc, coeff = 1.0):
    """Create operator for 4th integral of D P_n^{-1/2,-1/2}(2r^2-1) projected onto r P_n^{-1/2,1/2}(2r^2-1)."""
    if m != 1:
        raise RuntimeError("Operator only exists for m = 1!")

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,5)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*2.0/((2.0*n - 5.0)*(2.0*n - 3.0)*(2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*1.0/((n + 1.0)*(2.0*n - 3.0)*(2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -wb.worland_norm_row(n, m, -1)*3.0/(2.0*(n - 2.0)*(n + 1.0)*(2.0*n - 1.0))

    # Generate diagonal
    def d0(n):
        return -wb.worland_norm_row(n, m, 0)*3.0/(4.0*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*3.0*(2.0*n + 1.0)/(8.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*3.0*(2.0*n + 1.0)*(2.0*n + 3.0)/(16.0*n*(n + 1.0)**2*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -wb.worland_norm_row(n, m, 3)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)/(32.0*n*(n + 1.0)**2*(n + 2.0)**2*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -wb.worland_norm_row(n, m, 4)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)/(64.0*(n + 1.0)**2*(n + 2.0)**2*(n + 3.0)**2*(n + 4.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, m-1, bc)

def i4r_1d1(nr, m, bc, coeff = 1.0):
    """Create operator for 4th integral of 1/r D of r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*32.0*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -wb.worland_norm_row(n, m, -2)*96.0*m*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*24.0*(m + n - 1.0)*(4.0*m**2 - 4.0*m*n + 2.0*m - 4.0*n**2 + 4.0*n + 3.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate diagonal
    def d0(n):
        return -wb.worland_norm_row(n, m, 0)*16.0*m*(2.0*m**2 - 12.0*m*n - 12.0*n**2 + 7.0)/((m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -wb.worland_norm_row(n, m, 1)*6.0*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(4.0*m**2 - 4.0*m*n - 2.0*m - 4.0*n**2 - 4.0*n + 3.0)/((m + n)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -wb.worland_norm_row(n, m, 2)*6.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -wb.worland_norm_row(n, m, 3)*0.5*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, m, bc)

def i4laplh(nr, m, bc, coeff = 1.0):
    """Create operator for 4th integral of laplh r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*64.0*(m + n - 3.0)**2*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -wb.worland_norm_row(n, m, -2)*32.0*(m + n - 2.0)*(m + n - 1.0)*(4.0*m**2 - 9.0*m - 4.0*n**2 + 8.0*n + 5.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*16.0*(m + n - 1.0)*(4.0*m**3 - 20.0*m**2*n - 8.0*m**2 - 28.0*m*n**2 + 46.0*m*n + 28.0*m - 4.0*n**3 + 24.0*n**2 - 17.0*n - 15.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*16.0*(16.0*m**3*n + 6.0*m**3 - 36.0*m**2*n - 24.0*m**2 - 32.0*m*n**3 - 36.0*m*n**2 + 40.0*m*n + 21.0*m - 16.0*n**4 + 40.0*n**2 - 9.0)/((m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*4.0*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(24.0*m**2*n + 30.0*m**2 + 16.0*m*n**2 - 2.0*m*n - 45.0*m - 4.0*n**3 - 24.0*n**2 - 17.0*n + 15.0)/((m + n)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*2.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(8.0*m*n + 17.0*m + 4.0*n**2 + 8.0*n - 5.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*(n + 3.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, m, bc)

def i4lapl2h(nr, m, bc, coeff = 1.0):
    """Create operator for 4th integral of lapl2h r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*256.0*(m + n - 3.0)*(m + n - 2.0)**2*(m + n - 1.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*256.0*(m + n - 2.0)*(m + n - 1.0)*(4.0*m*n + m + 4.0*n**2 - 4.0*n - 3.0)/((m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*64.0*(24.0*m**2*n**2 + 36.0*m**2*n + 11.0*m**2 + 48.0*m*n**3 + 36.0*m*n**2 - 46.0*m*n - 27.0*m + 24.0*n**4 - 46.0*n**2 + 10.0)/((m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*64.0*(n + 2.0)*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(4.0*m*n + 5.0*m + 4.0*n**2 + 4.0*n - 3.0)/((m + n)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*16.0*(n + 2.0)*(n + 3.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 2)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 2)
    return radbc.constrain(mat, m, bc)

def i6(nr, m, bc, coeff = 1.0):
    """Create operator for 6th integral r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-6,7)
    nzrow = 5

    # Generate 6th subdiagonal
    def d_6(n):
        if m == 0:
            return wb.worland_norm_row(n, 0, -6)/((2.0*n - 11.0)*(2.0*n - 9.0)*(2.0*n - 7.0)*(2.0*n - 5.0)*(2.0*n - 3.0)*(2.0*n - 1.0))
        else:
            return wb.worland_norm_row(n, m, -6)*64.0*(m + n - 6.0)*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 12.0)*(m + 2.0*n - 11.0)*(m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -wb.worland_norm_row(n, m, -5)*384.0*m*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 11.0)*(m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return wb.worland_norm_row(n, m, -4)*96.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(10.0*m**2 - 4.0*m*n + 8.0*m - 4.0*n**2 + 16.0*n + 9.0)/((m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -wb.worland_norm_row(n, m, -3)*160.0*m*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(8.0*m**2 - 12.0*m*n + 18.0*m - 12.0*n**2 + 36.0*n + 37.0)/((m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*60.0*(m + n - 2.0)*(m + n - 1.0)*(16.0*m**4 - 64.0*m**3*n + 64.0*m**3 - 48.0*m**2*n**2 + 96.0*m**2*n + 236.0*m**2 + 32.0*m*n**3 - 96.0*m*n**2 - 40.0*m*n + 104.0*m + 16.0*n**4 - 64.0*n**3 - 40.0*n**2 + 208.0*n + 105.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -wb.worland_norm_row(n, m, -1)*48.0*m*(m + n - 1.0)*(8.0*m**4 - 80.0*m**3*n + 40.0*m**3 + 280.0*m**2 + 160.0*m*n**3 - 240.0*m*n**2 - 440.0*m*n + 260.0*m + 80.0*n**4 - 160.0*n**3 - 440.0*n**2 + 520.0*n + 537.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*4.0*(16.0*m**6 - 480.0*m**5*n + 960.0*m**4*n**2 + 1120.0*m**4 + 2560.0*m**3*n**3 - 7840.0*m**3*n + 480.0*m**2*n**4 - 5040.0*m**2*n**2 + 5614.0*m**2 - 960.0*m*n**5 + 5600.0*m*n**3 - 5180.0*m*n - 320.0*n**6 + 2800.0*n**4 - 5180.0*n**2 + 1125.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*12.0*m*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(8.0*m**4 - 80.0*m**3*n - 40.0*m**3 + 280.0*m**2 + 160.0*m*n**3 + 240.0*m*n**2 - 440.0*m*n - 260.0*m + 80.0*n**4 + 160.0*n**3 - 440.0*n**2 - 520.0*n + 537.0)/((m + n)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*15.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(16.0*m**4 - 64.0*m**3*n - 64.0*m**3 - 48.0*m**2*n**2 - 96.0*m**2*n + 236.0*m**2 + 32.0*m*n**3 + 96.0*m*n**2 - 40.0*m*n - 104.0*m + 16.0*n**4 + 64.0*n**3 - 40.0*n**2 - 208.0*n + 105.0)/(4.0*(m + n)*(m + n + 1.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*5.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(8.0*m**2 - 12.0*m*n - 18.0*m - 12.0*n**2 - 36.0*n + 37.0)/(2.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0))

    # Generate 4th superdiagonal
    def d4(n):
        return wb.worland_norm_row(n, m, 4)*3.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(10.0*m**2 - 4.0*m*n - 8.0*m - 4.0*n**2 - 16.0*n + 9.0)/(8.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0))

    # Generate 5rd superdiagonal
    def d5(n):
        return wb.worland_norm_row(n, m, 5)*3.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)/(8.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0)*(m + 2.0*n + 11.0))

    # Generate 6th superdiagonal
    def d6(n):
        return wb.worland_norm_row(n, m, 6)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*n + 11.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)*(2.0*m + 2.0*n + 11.0)/(64.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + n + 5.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0)*(m + 2.0*n + 11.0)*(m + 2.0*n + 12.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 3)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 3)
    return radbc.constrain(mat, m, bc)

def i6r_1d1(nr, m, bc, coeff = 1.0):
    """Create operator for 6th integral of 1/r D of r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-5,6)
    nzrow = 5

    # Generate 5th subdiagonal
    def d_5(n):
        return wb.worland_norm_row(n, m, -5)*128.0*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -wb.worland_norm_row(n, m, -4)*640.0*m*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*160.0*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(8.0*m**2 - 4.0*m*n + 6.0*m - 4.0*n**2 + 12.0*n + 7.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -wb.worland_norm_row(n, m, -2)*640.0*m*(m + n - 2.0)*(m + n - 1.0)*(2.0*m**2 - 4.0*m*n + 4.0*m - 4.0*n**2 + 8.0*n + 9.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*80.0*(m + n - 1.0)*(8.0*m**4 - 48.0*m**3*n + 24.0*m**3 - 32.0*m**2*n**2 + 32.0*m**2*n + 112.0*m**2 + 32.0*m*n**3 - 48.0*m*n**2 - 56.0*m*n + 36.0*m + 16.0*n**4 - 32.0*n**3 - 56.0*n**2 + 72.0*n + 45.0)/((m + 2*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate diagonal
    def d0(n):
        return -wb.worland_norm_row(n, m, 0)*16.0*m*(8.0*m**4 - 160.0*m**3*n + 80.0*m**2*n**2 + 260.0*m**2 + 480.0*m*n**3 - 920.0*m*n + 240.0*n**4 - 920.0*n**2 + 407.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -wb.worland_norm_row(n, m, 1)*20.0*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(8.0*m**4 - 48.0*m**3*n - 24.0*m**3 - 32.0*m**2*n**2 - 32.0*m**2*n + 112.0*m**2 + 32.0*m*n**3 + 48.0*m*n**2 - 56.0*m*n - 36.0*m + 16.0*n**4 + 32.0*n**3 - 56.0*n**2 - 72.0*n + 45.0)/((m + n)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -wb.worland_norm_row(n, m, 2)*40.0*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m**2 - 4.0*m*n - 4.0*m - 4.0*n**2 - 8.0*n + 9.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -wb.worland_norm_row(n, m, 3)*2.5*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(8.0*m**2 - 4.0*m*n - 6.0*m - 4.0*n**2 - 12.0*n + 7.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -wb.worland_norm_row(n, m, 4)*2.5*m*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0))

    # Generate 5rd superdiagonal
    def d5(n):
        return -wb.worland_norm_row(n, m, 5)*0.125*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 3)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 3)
    return radbc.constrain(mat, m, bc)

def i6laplh(nr, m, bc, coeff = 1.0):
    """Create operator for 6th integral of laplh r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-5,6)
    nzrow = 5

    # Generate 5th subdiagonal
    def d_5(n):
        return wb.worland_norm_row(n, m, -5)*256.0*(m + n - 5.0)**2*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 10.0)*(m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -wb.worland_norm_row(n, m, -4)*128.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(8.0*m**2 + 4.0*m*n - 33.0*m - 4.0*n**2 + 16.0*n + 9.0)/((m + 2.0*n - 9.0)*(m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*64.0*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(24.0*m**3 - 24.0*m**2*n - 64.0*m**2 - 60.0*m*n**2 + 230.0*m*n + 54.0*m - 12.0*n**3 + 104.0*n**2 - 183.0*n - 119.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -wb.worland_norm_row(n, m, -2)*128.0*(m + n - 2.0)*(m + n - 1.0)*(8.0*m**4 - 40.0*m**3*n - 10.0*m**3 - 48.0*m**2*n**2 + 196.0*m**2*n + 84.0*m**2 + 16.0*m*n**3 + 52.0*m*n**2 - 236.0*m*n - 157.0*m + 16.0*n**4 - 64.0*n**3 - 40.0*n**2 + 208.0*n + 105.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*32.0*(m + n - 1.0)*(8.0*m**5 - 152.0*m**4*n - 24.0*m**4 + 32.0*m**3*n**2 + 568.0*m**3*n + 388.0*m**3 + 448.0*m**2*n**3 - 272.0*m**2*n**2 - 1704.0*m**2*n - 636.0*m**2 + 272.0*m*n**4 - 944.0*m*n**3 - 832.0*m*n**2 + 2404.0*m*n + 1299.0*m + 16.0*n**5 - 240.0*n**4 + 360.0*n**3 + 800.0*n**2 - 891.0*n - 585.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*16.0*(96.0*m**5*n + 40.0*m**5 - 480.0*m**4*n**2 - 800.0*m**4*n - 360.0*m**4 - 960.0*m**3*n**3 + 400.0*m**3*n**2 + 3600.0*m**3*n + 1300.0*m**3 + 2400.0*m**2*n**3 + 1920.0*m**2*n**2 - 4600.0*m**2*n - 2880.0*m**2 + 576.0*m*n**5 + 1200.0*m*n**4 - 3360.0*m*n**3 - 4600.0*m*n**2 + 3108.0*m*n + 2035.0*m + 192.0*n**6 - 1680.0*n**4 + 3108.0*n**2 - 675.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*8.0*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(120.0*m**4*n + 160.0*m**4 - 160.0*m**3*n**2 - 760.0*m**3*n - 900.0*m**3 - 480.0*m**2*n**3 - 1120.0*m**2*n**2 + 1040.0*m**2*n + 2240.0*m**2 - 192.0*m*n**4 + 16.0*m*n**3 + 1912.0*m*n**2 + 804.0*m*n - 2190.0*m + 16.0*n**5 + 240.0*n**4 + 360.0*n**3 - 800.0*n**2 - 891.0*n + 585.0)/((m + n)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*8.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(40.0*m**3*n + 90.0*m**3 - 100.0*m**2*n - 280.0*m**2 - 48.0*m*n**3 - 244.0*m*n**2 - 156.0*m*n + 365.0*m - 16.0*n**4 - 64.0*n**3 + 40.0*n**2 + 208.0*n - 105.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(60.0*m**2*n + 190.0*m**2 + 24.0*m*n**2 + 22.0*m*n - 237.0*m - 12.0*n**3 - 104.0*n**2 - 183.0*n + 119.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    # Generate 4th superdiagonal
    def d4(n):
        return wb.worland_norm_row(n, m, 4)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(12.0*m*n + 49.0*m + 4.0*n**2 + 16.0*n - 9.0)/(2.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0))

    # Generate 5rd superdiagonal
    def d5(n):
        return wb.worland_norm_row(n, m, 5)*(n + 5.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*n + 9.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)*(2.0*m + 2.0*n + 9.0)/(4.0*(m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + n + 4.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0)*(m + 2.0*n + 9.0)*(m + 2.0*n + 10.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 3)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 3)
    return radbc.constrain(mat, m, bc)

def i6lapl2h(nr, m, bc, coeff = 1.0):
    """Create operator for 6th integral of lapl2h r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-4,5)
    nzrow = 5

    # Generate 4th subdiagonal
    def d_4(n):
        return wb.worland_norm_row(n, m, -4)*1024.0*(m + n - 5.0)*(m + n - 4.0)**2*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 8.0)*(m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -wb.worland_norm_row(n, m, -3)*1024.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(2.0*m**2 - 2.0*m*n - 7.0*m - 4.0*n**2 + 12.0*n + 7.0)/((m + 2.0*n - 7.0)*(m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*256.0*(m + n - 2.0)*(m + n - 1.0)*(4.0*m**4 - 32.0*m**3*n - 28.0*m**3 - 60.0*m**2*n**2 + 240.0*m**2*n + 125.0*m**2 - 8.0*m*n**3 + 244.0*m*n**2 - 530.0*m*n - 321.0*m + 16.0*n**4 - 24.0*n**3 - 236.0*n**2 + 390.0*n + 250.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*256.0*(m + n - 1.0)*(24.0*m**4*n + 8.0*m**4 - 24.0*m**3*n**2 - 196.0*m**3*n - 92.0*m**3 - 136.0*m**2*n**3 + 24.0*m**2*n**2 + 706.0*m**2*n + 358.0*m**2 - 104.0*m*n**4 + 388.0*m*n**3 + 362.0*m*n**2 - 1101.0*m*n - 562.0*m - 16.0*n**5 + 160.0*n**4 - 200.0*n**3 - 520.0*n**2 + 531.0*n + 360.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*64.0*(240.0*m**4*n**2 + 400.0*m**4*n + 136.0*m**4 + 320.0*m**3*n**3 - 800.0*m**3*n**2 - 2432.0*m**3*n - 1000.0*m**3 - 240.0*m**2*n**4 - 2400.0*m**2*n**3 - 816.0*m**2*n**2 + 5000.0*m**2*n + 2375.0*m**2 - 480.0*m*n**5 - 1200.0*m*n**4 + 3232.0*m*n**3 + 5000.0*m*n**2 - 3130.0*m*n - 2375.0*m - 160.0*n**6 + 1616.0*n**4 - 3130.0*n**2 + 684.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*64.0*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(80.0*m**3*n**2 + 280.0*m**3*n + 236.0*m**3 + 120.0*m**2*n**3 + 180.0*m**2*n**2 - 618.0*m**2*n - 939.0*m**2 + 24.0*m*n**4 - 252.0*m*n**3 - 962.0*m*n**2 - 61.0*m*n + 1093.0*m - 16.0*n**5 - 160.0*n**4 - 200.0*n**3 + 520.0*n**2 + 531.0*n - 360.0)/((m + n)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*16.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(60.0*m**2*n**2 + 320.0*m**2*n + 419.0*m**2 + 72.0*m*n**3 + 316.0*m*n**2 + 58.0*m*n - 711.0*m + 16.0*n**4 + 24.0*n**3 - 236.0*n**2 - 390.0*n + 250.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*16.0*(n + 4.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(6.0*m*n + 19.0*m + 4.0*n**2 + 12.0*n - 7.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0))

    # Generate 4th superdiagonal
    def d4(n):
        return wb.worland_norm_row(n, m, 4)*4.0*(n + 4.0)*(n + 5.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)*(2.0*m + 2.0*n + 7.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + n + 3.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0)*(m + 2.0*n + 7.0)*(m + 2.0*n + 8.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 3)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 3)
    return radbc.constrain(mat, m, bc)

def i6lapl3h(nr, m, bc, coeff = 1.0):
    """Create operator for 6th integral of lapl3h r^m P_n^{-1/2,m-1/2}(2r^2-1)."""

    ns = np.arange(0, nr+3)
    offsets = np.arange(-3,4)
    nzrow = 5

    # Generate 3rd subdiagonal
    def d_3(n):
        return wb.worland_norm_row(n, m, -3)*4096.0*(m + n - 5.0)*(m + n - 4.0)*(m + n - 3.0)**2*(m + n - 2.0)*(m + n - 1.0)/((m + 2.0*n - 6.0)*(m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return wb.worland_norm_row(n, m, -2)*6144.0*(m + n - 4.0)*(m + n - 3.0)*(m + n - 2.0)*(m + n - 1.0)*(4.0*m*n + m + 4.0*n**2 - 8.0*n - 5.0)/((m + 2.0*n - 5.0)*(m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row(n, m, -1)*3072.0*(m + n - 3.0)*(m + n - 1.0)*(20.0*m**2*n**2 + 30.0*m**2*n + 9.0*m**2 + 40.0*m*n**3 - 10.0*m*n**2 - 109.0*m*n - 48.0*m + 20.0*n**4 - 40.0*n**3 - 59.0*n**2 + 79.0*n + 48.0)/((m + 2.0*n - 4.0)*(m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0))

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row(n, m, 0)*512.0*(160.0*m**3*n**3 + 600.0*m**3*n**2 + 656.0*m**3*n + 195.0*m**3 + 480.0*m**2*n**4 + 1200.0*m**2*n**3 - 480.0*m**2*n**2 - 2370.0*m**2*n - 954.0*m**2 + 480.0*m*n**5 + 600.0*m*n**4 - 2272.0*m*n**3 - 2370.0*m*n**2 + 1930.0*m*n + 1245.0*m + 160.0*n**6 - 1136.0*n**4 + 1930.0*n**2 - 414.0)/((m + 2.0*n - 3.0)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return wb.worland_norm_row(n, m, 1)*768.0*(n + 3.0)*(2.0*n + 1.0)*(2.0*m + 2.0*n + 1.0)*(20.0*m**2*n**2 + 70.0*m**2*n + 59.0*m**2 + 40.0*m*n**3 + 110.0*m*n**2 - 9.0*m*n - 127.0*m + 20.0*n**4 + 40.0*n**3 - 59.0*n**2 - 79.0*n + 48.0)/((m + n)*(m + 2.0*n - 2.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return wb.worland_norm_row(n, m, 2)*384.0*(n + 3.0)*(n + 4.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(4.0*m*n + 9.0*m + 4.0*n**2 + 8.0*n - 5.0)/((m + n)*(m + n + 1.0)*(m + 2.0*n - 1.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return wb.worland_norm_row(n, m, 3)*64.0*(n + 3.0)*(n + 4.0)*(n + 5.0)*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*m + 2.0*n + 1.0)*(2.0*m + 2.0*n + 3.0)*(2.0*m + 2.0*n + 5.0)/((m + n)*(m + n + 1.0)*(m + n + 2.0)*(m + 2.0*n + 1.0)*(m + 2.0*n + 2.0)*(m + 2.0*n + 3.0)*(m + 2.0*n + 4.0)*(m + 2.0*n + 5.0)*(m + 2.0*n + 6.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 3)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 3)
    return radbc.constrain(mat, m, bc)

def divr(nr, m, bc, coeff = 1.0):
    """Create operator for 1/r r^m P_n^{-1/2,m-3/2}(2r^2-1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,1)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        return wb.worland_norm_row_l_1(n,m,-1)*n/(2.0*n + m - 2.0)

    # Generate diagonal
    def d0(n):
        return wb.worland_norm_row_l_1(n,m,0)*(2.0*n + 2.0*m - 1.0)/(2.0*(2.0*n + m))

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    mat = radbc.restrict_eye(mat.shape[0], 'rt', 1)*mat*radbc.restrict_eye(mat.shape[1], 'cr', 1)
    return radbc.constrain(mat, m, bc)

def qid(nr, m, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = spsp.coo_matrix((nr,nr))
    if coeff != 1.0:
        mat.data = coeff*np.ones((nr-q))
    else:
        mat.data = np.ones((nr-q))
    mat.row = np.arange(q,nr)
    mat.col = mat.row
    return radbc.constrain(mat, m, bc)

def integral(nr, m):
    """Compute definite integral of the expansion"""

    print("IMPLEMENTATION OF INTEGRAL IS WRONG!")
    mat = spsp.lil_matrix((1,nr))
    mat[0,::2] = [4.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,nr,2)]
    mat[0,0] = mat[0,0]/2.0

    return mat.tocoo()

def stencil(nr, m, bc, make_square):
    """Create a galerkin stencil matrix"""

    mat = qid(nr, m, 0, radbc.no_bc())

    if make_square:
        bc['rb'] = bc['rt']
    bc['rt'] = 0

    return radbc.constrain(mat, m, bc)