"""Module provides functions to generate sparse operators in a cartesian box with two periodic directions."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc


def zblk(nx, bc):
    """Create a block of zeros"""

    mat = spsp.lil_matrix((nx,nx))
    return c1dbc.constrain(mat,bc)

def d1(nx, bc, coeff = 1.0, cscale = 1.0, zr = 1):
    """Create operator for 1st derivative"""

    row = [2*j for j in range(0,nx)]
    mat = spsp.lil_matrix((nx,nx))
    for i in range(0,nx-1):
        mat[i,i+1:nx:2] = row[i+1:nx:2]
    mat[-zr:,:] = 0

    mat = coeff*cscale*mat
    return c1dbc.constrain(mat, bc, location = 'b')

def d2(nx, bc, coeff = 1.0, cscale = 1.0, zr = 2):
    """Create operator for 2nd derivative"""

    mat = spsp.lil_matrix((nx,nx))
    for i in range(0,nx-2):
        mat[i,i+2:nx:2] = [j*(j**2 - i**2) for j in range(0,nx)][i+2:nx:2]
    mat[-zr:,:] = 0

    mat = coeff*cscale**2*mat
    return c1dbc.constrain(mat, bc, location = 'b')

def d4(nx, bc, coeff = 1.0, cscale = 1.0, zr = 4):
    """Create operator for 4th derivative"""

    mat_d2 = d2(nx + 4, c1dbc.no_bc())
    mat = mat_d2*mat_d2
    mat = mat[0:-4, 0:-4]
    mat[-zr:,:] = 0
    
    mat = coeff*cscale**4*mat
    return c1dbc.constrain(mat, bc, location = 'b')

def laplh(nx, k, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for horizontal laplacian"""

    mat = d2(nx, bc, cscale = cscale) - k**2*sid(nx,2,bc)

    mat = coeff*mat
    return c1dbc.constrain(mat, bc, location = 'b')

def lapl2h(nx, k, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for horizontal bilaplacian"""

    mat = d4(nx, bc, cscale = cscale) - 2.0*k**2*sid(nx,4,bc)*d2(nx,bc, cscale = cscale) + k**4*sid(nx,4,bc) 

    mat = coeff*mat
    return c1dbc.constrain(mat, bc, location = 'b')

def i1d1(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create a quasi identity block of order 1"""

    return qid(nx,1, bc, coeff*cscale)

def i2d2(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create a quasi identity block of order 2"""

    return qid(nx,2, bc, coeff*cscale**2)

def i1(nx, bc, coeff = 1.0):
    """Create operator for 1st integral in x of T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-1,2,2)
    nzrow = 0

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(2.0*n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i2(nx, bc, coeff = 1.0):
    """Create operator for 2nd integral in x of T_n(x)"""
    
    ns = np.arange(0, nx, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -1.0/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 1.0/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i2x1(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 2nd integral in x of x T_n(x)"""
    
    ns = np.arange(0, nx, 1)
    offsets = np.arange(-3,4,2)
    nzrow = 1

    # Generate 3rd subdiagonal
    def d_3(n):
        return 1.0/(8.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -1.0/(8.0*n*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(8.0*n*(n - 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return 1.0/(8.0*n*(n + 1.0))

    ds = [d_3, d_1, d1, d3]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*cscale*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i2d1(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 2nd integral in x of T_n(x)"""
    
    ns = np.arange(0, nx, 1)
    offsets = np.arange(-1,2,2)
    nzrow = 1

    # Generate 1st subdiagonal
    def d_1(n):
        return 1.0/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -1.0/(2.0*n)

    ds = [d_1, d1]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*cscale*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i2lapl(nx, k, l, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 2nd integral in x of Laplacian T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(k**2 + l**2)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (k**2 + l**2 + (2.0*n**2 - 2.0)*cscale**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(k**2 + l**2)/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i2laplh(nx, k, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 2nd integral in x of horizontal Laplacian T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    # Generate 2nd subdiagonal
    def d_2(n):
        return -k**2/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (k**2 + (2.0*n**2 - 2.0)*cscale**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -k**2/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4(nx, bc, coeff = 1.0):
    """Create operator for 4th integral in x of T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3
   
    # Generate 4th subdiagonal
    def d_4(n):
        return 1.0/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -1.0/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

    # Generate main diagonal
    def d0(n):
        return 3.0/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -1.0/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0)) 

    # Generate 4th superdiagonal
    def d4(n):
        return 1.0/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4d1(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of D_x T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-3,4,2)
    nzrow = 3

    # Generate 3rd subdiagonal
    def d_2(n):
        return 1.0/(8.0*n*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0/(8.0*n*(n - 2.0)*(n + 1.0)) 

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0/(8.0*n*(n - 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -1.0/(8.0*n*(n + 1.0)*(n + 2.0)) 

    ds = [d_2, d_1, d1, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*cscale*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4d2(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of D_x^2 T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 3

    # Generate 2nd subdiagonal
    def d_2(n):
        return 1.0/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -1.0/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return 1.0/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*cscale**2*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4d4(nx, bc, coeff = 1.0, cscale = 1.0):
    """Create a quasi identity block of order 4"""

    return qid(nx,4, bc, coeff*cscale**4)

def i4lapl(nx, k, l, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of Laplacian T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3
   
    # Generate 4th subdiagonal
    def d_4(n):
        return -(k**2 + l**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (k**2 + l**2 + (n**2 - 2.0*n - 3.0)*cscale**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

    # Generate main diagonal
    def d0(n):
        return -(3.0*k**2 + 3.0*l**2 + (4.0*n**2 - 16.0)*cscale**2)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (k**2 + l**2 + (n**2 + 2.0*n - 3.0)*cscale**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(k**2 + l**2)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4laplh(nx, k, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of horizontal Laplacian T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3
   
    # Generate 4th subdiagonal
    def d_4(n):
        return -k**2/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (k**2 + (n**2 - 2.0*n - 3.0)*cscale**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

    # Generate main diagonal
    def d0(n):
        return -(3.0*k**2 + (4.0*n**2 - 16.0)*cscale**2)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (k**2 + (n**2 + 2.0*n - 3.0)*cscale**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -k**2/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4lapl2(nx, k, l, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of Laplacian^2 T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return (k**2 + l**2)**2/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(k**2 + l**2)*(k**2 + l**2 + (2.0*n**2 - 4.0*n - 6.0)*cscale**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

    # Generate main diagonal
    def d0(n):
        return (3.0*k**4 + 6.0*k**2*l**2 + 3.0*l**4 + (8.0*k**2*n**2 - 32.0*k**2 + 8.0*l**2*n**2 - 32.0*l**2)*cscale**2 + (8.0*n**4 - 40.0*n**2 + 32.0)*cscale**4)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(k**2 + l**2)*(k**2 + l**2 + (2.0*n**2 + 4.0*n - 6.0)*cscale**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0)) 

    # Generate 4th superdiagonal
    def d4(n):
        return (k**2 + l**2)**2/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def i4lapl2h(nx, k, bc, coeff = 1.0, cscale = 1.0):
    """Create operator for 4th integral in x of horizontal Laplacian^2 T_n(x)"""

    ns = np.arange(0, nx, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    # Generate 4th subdiagonal
    def d_4(n):
        return k**4/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -k**2*(k**2 + (2.0*n**2 - 4.0*n - 6.0)*cscale**2)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)) 

    # Generate main diagonal
    def d0(n):
        return (3.0*k**4 + (8.0*k**2*n**2 - 32.0*k**2)*cscale**2 + (8.0*n**4 - 40.0*n**2 + 32.0)*cscale**4)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -k**2*(k**2 + (2.0*n**2 + 4.0*n - 6.0)*cscale**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0)) 

    # Generate 4th superdiagonal
    def d4(n):
        return k**4/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags = utils.build_diagonals(ns, nzrow, ds, offsets)

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def qid(nx, q, bc, coeff = 1.0):
    """Create a quasi identity block of order q"""

    offsets = [0]
    diags = [[0]*q + [1]*(nx-q)]

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc)

def sid(nx, s, bc, coeff = 1.0):
    """Create a identity block with last s rows zeroed"""

    offsets = [0]
    diags = [[1]*(nx-s) + [0]*s]

    mat = coeff*spsp.diags(diags, offsets)
    return c1dbc.constrain(mat, bc, location = 'b')

def stencil(nx, bc):
    """Create a galerkin stencil matrix"""

    bc['rt'] = 0
    return c1dbc.stencil(nx, bc)

def avg(nx):
    """Compute the average of the expansion"""

    mat = spsp.lil_matrix((1,nx))
    mat[0,::2] = [2.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,nx,2)]
    mat[0,0] = mat[0,0]/2.0

    return mat

def integral(nr, cscale = 1.0):
    """Compute the definite integral of the expansion"""

    mat = spsp.lil_matrix((1,nr))
    mat[0,::2] = [4.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,nr,2)]
    mat[0,0] = mat[0,0]/2.0

    return mat

def surfaceFlux(nx, cscale = 1.0):
    """Compute the flux through boundary"""

    mat = c1dbc.constrain(spsp.lil_matrix((1, nx)), {0:12})
    mat = 2.0*cscale*mat

    return mat
