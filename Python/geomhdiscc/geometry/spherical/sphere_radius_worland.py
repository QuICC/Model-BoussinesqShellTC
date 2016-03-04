"""Module provides functions to generate sparse operators for the radial direction in a sphere using r^ P_n^{-1/2,l - 1/2}(2r^2 -1)."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.sphere_radius_worland_boundary as radbc


def zblk(nr, l, bc):
    """Create a block of zeros"""

    mat = spsp.coo_matrix((nr,nr))
    return radbc.constrain(mat,l,bc)

def i1(nr, l, bc, coeff = 1.0):
    """Create operator for 1st integral r^l P_n^{-1/2,l-1/2}(2r^2-1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-1,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 2.0*(l + n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -2.0*l/((l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)/(2.0*(l + n)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Remove extra column and unused boundary row
    bc['cr'] = bc.get('cr',0) + 1
    bc['rt'] = bc.get('rt',0) + 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 4.0*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -8.0*l*(l + n - 1.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return 2.0*(2.0*l**2 - 4.0*l*n - 4.0*n**2 + 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 2.0*l*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/(4.0*(l + n)*(l + n + 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, has_wrap = False)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-4,5)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return 16.0*(l + n - 4.0)*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 8.0)*(l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -64.0*l*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)/((l + 2.0*n - 7.0)*(l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return 16.0*(l + n - 2.0)*(l + n - 1.0)*(6.0*l**2 - 4.0*l*n + 4.0*l - 4.0*n**2 + 8.0*n + 5.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -16.0*l*(l + n - 1.0)*(4.0*l**2 - 12.0*l*n + 6.0*l - 12.0*n**2 + 12.0*n + 17.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate diagonal
    def d0(n):
        return 2.0*(8.0*l**4 - 96.0*l**3*n - 48.0*l**2*n**2 + 100.0*l**2 + 96.0*l*n**3 - 120.0*l*n + 48.0*n**4 - 120.0*n**2 + 27.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 4.0*l*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)*(4.0*l**2 - 12.0*l*n - 6.0*l - 12.0*n**2 - 12.0*n + 17.0)/((l + n)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(6.0*l**2 - 4.0*l*n - 4.0*l - 4.0*n**2 - 8.0*n + 5.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return l*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/((l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0))

    # Generate 4th superdiagonal
    def d4(n):
        return (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)*(2.0*n + 7.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)*(2.0*l + 2.0*n + 7.0)/(16.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + n + 3.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0)*(l + 2.0*n + 7.0)*(l + 2.0*n + 8.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i2lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 2nd integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr+1)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Remove extra column and unused boundary row
    bc['cr'] = bc.get('cr',0) + 1
    bc['rt'] = bc.get('rt',0) + 1

    # Generate 1st subdiagonal
    def d_1(n):
        return 8.0*(l + n - 1.0)*(2.0*l + 2.0*n - 1.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate main diagonal
    def d0(n):
        return 8.0*(4.0*l*n + 4.0*n**2 - 1.0)/((l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 2.0*(2.0*n + 1.0)**2*(2.0*l + 2.0*n + 1.0)/((l + n)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4lapl(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral of Laplacian r^l P_n^{-1/2,l-1/2}(2r^2 -1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return 32.0*(l + n - 3.0)*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)/((l + 2.0*n - 6.0)*(l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -32.0*(l + n - 2.0)*(l + n - 1.0)*(4.0*l**2 - 6.0*l - 4.0*n**2 + 8.0*n + 5.0)/((l + 2.0*n - 5.0)*(l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 8.0*(l + n - 1.0)*(8.0*l**3 - 40.0*l**2*n - 4.0*l**2 - 56.0*l*n**2 + 80.0*l*n + 62.0*l - 8.0*n**3 + 36.0*n**2 - 22.0*n - 21.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate main diagonal
    def d0(n):
        return 16.0*(16.0*l**3*n + 4.0*l**3 - 24.0*l**2*n - 24.0*l**2 - 32.0*l*n**3 - 24.0*l*n**2 + 40.0*l*n + 14.0*l - 16.0*n**4 + 40.0*n**2 - 9.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 2.0*(2.0*n + 1.0)*(2.0*l + 2.0*n + 1.0)*(48.0*l**2*n + 48.0*l**2 + 32.0*l*n**2 + 8.0*l*n - 84.0*l - 8.0*n**3 - 36.0*n**2 - 22.0*n + 21.0)/((l + n)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 2.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(8.0*l*n + 14.0*l + 4.0*n**2 + 8.0*n - 5.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return (2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*n + 5.0)**2*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)*(2.0*l + 2.0*n + 5.0)/(2.0*(l + n)*(l + n + 1.0)*(l + n + 2.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0)*(l + 2.0*n + 5.0)*(l + 2.0*n + 6.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def i4lapl2(nr, l, bc, coeff = 1.0):
    """Create operator for 4th integral bilaplacian r^l P_n^{-1/2, l-1/2}(2r^2 - 1)."""

    ns = np.arange(0, nr)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return 64.0*(l + n - 2.0)*(l + n - 1.0)*(2.0*l + 2.0*n - 5.0)*(2.0*l + 2.0*n - 3.0)/((l + 2.0*n - 4.0)*(l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 128.0*(l + n - 1.0)*(2.0*l + 2.0*n - 3.0)*(4.0*l*n + 4.0*n**2 - 4.0*n - 3.0)/((l + 2.0*n - 3.0)*(l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return 96.0*(2.0*n + 1.0)*(2.0*l + 2.0*n - 1.0)*(4.0*l*n + 2.0*l + 4.0*n**2 - 5.0)/((l + 2.0*n - 2.0)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 32.0*(2.0*n + 1.0)*(2.0*n + 3.0)*(2.0*l + 2.0*n + 1.0)*(4.0*l*n + 4.0*l + 4.0*n**2 + 4.0*n - 3.0)/((l + n)*(l + 2.0*n - 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 4.0*(2.0*n + 1.0)*(2.0*n + 3.0)**2*(2.0*n + 5.0)*(2.0*l + 2.0*n + 1.0)*(2.0*l + 2.0*n + 3.0)/((l + n)*(l + n + 1.0)*(l + 2.0*n + 1.0)*(l + 2.0*n + 2.0)*(l + 2.0*n + 3.0)*(l + 2.0*n + 4.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets, format = 'coo')
    return radbc.constrain(mat, l, bc)

def qid(nr, l, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""
    
    mat = spsp.coo_matrix((nr,nr))
    if coeff != 1.0:
        mat.data = coeff*np.ones((nr-q))
    else:
        mat.data = np.ones((nr-q))
    mat.row = np.arange(q,nr)
    mat.col = mat.row
    return radbc.constrain(mat, l, bc)

def stencil(nr, l, bc, make_square):
    """Create a galerkin stencil matrix"""

    mat = qid(nr, l, 0, radbc.no_bc())

    if not make_square:
        bc['rt'] = 0

    return radbc.constrain(mat, l, bc)
