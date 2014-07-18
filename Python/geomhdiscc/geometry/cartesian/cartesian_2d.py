"""Module provides functions to generate sparse operators in a cartesian box with a single periodic dimension."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_boundary_2d as c2dbc


def convert_bc(bc):
    """Convert boundary dictionary into x and z kronecker product boundaries"""

    if bc['x'][0] < 0:
        bcx = bc['x']
    else:
        bcx = c1d.c1dbc.no_bc()
        bcx['r'] = bc['x'].get('r',0)

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        bcz['r'] = bc['z'].get('r',0)

    return (bcx, bcz)

def zblk(nx, nz, bc):
    """Create a block of zeros"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,bcz),c1d.zblk(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 0, 0, bc)

def i2j1d0d1(nx, nz, bc, coeff = 1.0):
    """Create a quasi identity block of order 2,2"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,1,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 1, bc)

def i2j2d0d2(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 2nd integrazl in z of 2nd derivative"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2d2d2(nx, nz, bc, coeff = 1.0):
    """Create a quasi identity block of order 2,2"""

    return qid(nx,nz,2,2, bc, coeff)

def i2j0(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,0,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 0, bc)

def i2j1(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 1st integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 1, bc)

def i2j2(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j0laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 1st integral in z of of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,0,bcz), c1d.i2laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 0, bc)

def i2j1laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 1st integral in z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i2laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 1, bc)

def i2j2laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), c1d.i2laplh(nx,k,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2lapl(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), c1d.i2laplh(nx,k,bcx)) + spsp.kron(c1d.i2d2(nz,bcz), c1d.i2(nx,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i4j1(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 1st integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j2(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 2nd integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 2, bc)

def i4j1d0d1(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 1st integral of 1st derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,1,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j2d0d2(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x, 2nd integral of 2nd derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 2, bc)

def i4j4(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def i4j1laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x 1st integral of z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i4laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j2laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x 2nd integral of z of horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i4laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 2, bc)

def i4j4lapl(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4laplh(nx,k,bcx)) + spsp.kron(c1d.i4d2(nz,bcz), c1d.i4(nx,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def i4j1lapl2h(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 1st integral of z of horizontal Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i4lapl2h(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j4lapl2(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4lapl2h(nx,k,bcx)) + 2*spsp.kron(c1d.i4d2(nz,bcz), c1d.i4laplh(nx,k,bcx)) + spsp.kron(c1d.i4d4(nz,bcz), c1d.i4(nx,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def qid(nx, nz, qx, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qx in x in qz in z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), c1d.qid(nx,qx,bcx))
    return c2dbc.constrain(mat, nx, nz, qx, qz, bc)
