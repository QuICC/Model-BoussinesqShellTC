"""Module provides functions to generate sparse operators for the radial direction in a cylinder."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.cylinder_radius_boundary as radbc


def zblk(nr, parity, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nr,nr))
    return radbc.constrain(mat,parity,bc)

def x1(nr, parity, bc, coeff = 1.0, zr = 0):
    """Create operator for x multiplication"""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = -1

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/2.0

    # Generate 1st superdiagonal
    def d1(n):
        return 1.0/2.0

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    if zr > 0:
        mat = mat.tolil()
        mat[-zr:,:] = 0
        mat = mat.tocsr()
    return radbc.constrain(mat, parity, bc)

def x2(nr, parity, bc, coeff = 1.0, zr = 0):
    """Create operator for x^2 multiplication"""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = -1

    # Generate 2nd subdiagonal
    def d_1(n):
        return 1.0/4.0

    # Generate main diagonal
    def d0(n):
        return 1.0/2.0

    # Generate 2nd superdiagonal
    def d1(n):
        return 1.0/4.0

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def d1(nr, parity, bc, coeff = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    row = [2*j for j in range(parity,2*nr,2)]
    mat = spsp.lil_matrix((nr,nr))
    for i in range(0,nr):
        mat[i,i+(parity+1)%2:] = row[i+(parity+1)%2:]
    if zr > 0:
        mat[-zr:,:] = 0

    mat = coeff*mat
    return radbc.constrain(mat, parity, bc)

def x1d1(nr, parity, bc, coeff = 1.0, zr = 1):
    """Create operator for x times derivative"""

    mat = x1(nr, (parity+1)%2, radbc.no_bc(), coeff, zr = zr)*d1(nr, parity, radbc.no_bc(), coeff, zr = zr)
    return radbc.constrain(mat, parity, bc)

def x1div(nr, parity, bc, coeff = 1.0, zr = 1):
    """Create operator for x times radial divergence"""

    mat = sid(nr, parity, zr, radbc.no_bc(), coeff) + x1(nr, (parity+1)%2, radbc.no_bc(), coeff, zr = zr)*d1(nr, parity, radbc.no_bc(), coeff, zr = zr)
    return radbc.constrain(mat, parity, bc)

def i1(nr, parity, bc, coeff = 1.0):
    """Create operator for 1st integral T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i1x1d1(nr, parity, bc, coeff = 1.0):
    """Create operator for 1st integral x 1st derivative T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return (n - 1.0)/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return (n + 1.0)/(2.0*n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i1x1div(nr, parity, bc, coeff = 1.0):
    """Create operator for 1st integral x radial divergence T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(0,2)
    else:
        offsets = np.arange(-1,1)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/2.0

    # Generate 1st superdiagonal
    def d1(n):
        return 1.0/2.0

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i1x1(nr, parity, bc, coeff = 1.0):
    """Create operator for 1st integral of x T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(4.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i1x2(nr, parity, bc, coeff = 1.0):
    """Create operator for 1st integral of x^2 T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 0

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(8.0*n)

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(8.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(8.0*n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -1.0/(8.0*n)

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_1(n):
        return 1.0/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -1.0/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d1(n):
        return 1.0/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x1(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -1.0/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 1.0/(8.0*n*(n + 1.0))

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x2d2(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 2nd derivative T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return (n - 3.0)*(n - 2.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n): 
        return (n**2 - 3.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (n + 2.0)*(n + 3.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x2d1(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 1st derivative T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
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
        return -(n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n + 3.0)/(8.0*n*(n + 1.0))

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x2(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 T_n(x)."""

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
    return radbc.constrain(mat, parity, bc)

def i2x2div(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 radial divergence T_n(x)."""

    ns = np.arange((parity+1)%2, 2*nr, 2)
    if parity == 0:
        offsets = np.arange(-1,3)
    else:
        offsets = np.arange(-2,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return (n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (n + 2.0)/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n + 2.0)/(8.0*n*(n + 1.0))

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets, (parity+1)%2)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x2laplh(nr, m, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^2 Laplacian T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return  -(m - n + 2.0)*(m + n - 2.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (m**2 + n**2 - 2.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(m - n - 2.0)*(m + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x3laplhx_1(nr, m, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^3 Laplacian 1/x T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return  -(m**2 - n**2 + 6.0*n - 8.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (m**2 + n**2 - 4.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(m**2 - n**2 - 6.0*n - 8.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x3d1(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^3 first derivative T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return (n - 4.0)/(16.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (n - 2.0)*(n + 2.0)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return 3.0/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -(n + 4.0)/(16.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i2x3d1x_2(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral of x^3 first derivative 1/x^2 T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-1,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return  (n - 4.0)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return 3.0/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(n + 4.0)/(4.0*n*(n + 1.0))

    ds = [d_1, d0, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i4(nr, parity, bc, coeff = 1.0):
    """Create operator for 2nd integral T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -1.0/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return 3.0/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 1.0/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i4x4(nr, parity, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 T_n(x)."""

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
    return radbc.constrain(mat, parity, bc)

def i4x4laplh(nr, m, parity, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-3,4)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(m - n + 6.0)*(m + n - 6.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (n - 5.0)*(m**2 + n**2 + 2.0*n - 24.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (m**2*n + 17.0*m**2 - n**3 + 27.0*n**2 - 8.0*n - 372.0)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -(m**2*n**2 - 19.0*m**2 + n**4 - 43.0*n**2 + 396.0)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (m**2*n - 17.0*m**2 - n**3 - 27.0*n**2 - 8.0*n + 372.0)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (n + 5.0)*(m**2 + n**2 - 2.0*n - 24.0)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(m - n - 6.0)*(m + n + 6.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def i4x4lapl2h(nr, m, parity, bc, coeff = 1.0):
    """Create operator for 4th integral of x^4 Laplacian^2 T_n(x)."""

    ns = np.arange(parity, 2*nr, 2)
    offsets = np.arange(-2,3)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return (m - n + 4.0)*(m - n + 6.0)*(m + n - 6.0)*(m + n - 4.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(m - n + 4.0)*(m + n - 4.0)*(m**2 + n**2 - 2.0*n - 12.0)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*m**4 + 2.0*m**2*n**2 - 68.0*m**2 + 3.0*n**4 - 68.0*n**2 + 416.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(m - n - 4.0)*(m + n + 4.0)*(m**2 + n**2 + 2.0*n - 12.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (m - n - 6.0)*(m - n - 4.0)*(m + n + 4.0)*(m + n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def qid(nr, parity, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    offsets = [0]
    diags = [[0]*q + [1]*(nr-q)]

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)

def sid(nr, parity, s, bc, coeff = 1.0):
    """Create a identity block with last s rows zeroed"""

    offsets = [0]
    diags = [[1]*(nr-s) + [0]*s]

    mat = coeff*spsp.diags(diags, offsets)
    return radbc.constrain(mat, parity, bc)
