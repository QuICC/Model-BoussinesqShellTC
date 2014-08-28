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
        for key, val in bc['x'].items():
            if key != 0:
                bcx[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcx, bcz)

def d1(nx, nz, bc, coeff = 1.0, sx = 1, sz = 0):
    """Create operator for the 1st Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz, sz, bcz), c1d.d1(nx, bcx, zr = sx))
    return c2dbc.constrain(mat, nx, nz, sx, sz, bc, location = 'b')

def e1(nx, nz, bc, coeff = 1.0, zscale = 1.0, sx = 0, sz = 1):
    """Create operator for the 1st Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.d1(nz, bcz, zscale, zr = sz), c1d.sid(nx, sx, bcx))
    return c2dbc.constrain(mat, nx, nz, sx, sz, bc, location = 'b')

def lapl(nx, nz, k, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for the 2nd X derivative and 2nd Z derivative T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.sid(nz, 2, bcz), c1d.laplh(nx, k, bcx))
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.d2(nz, bcz, zscale**2), c1d.sid(nx, 2, bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc, location = 'b')

def laplh(nx, nz, k, sz, bc, coeff = 1.0):
    """Create operator for the horizontal Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz,sz,bcz), c1d.laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, sz, bc, location = 'b')

def lapl2h(nx, nz, k, sz, bc, coeff = 1.0):
    """Create operator for the horizontal bilaplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz,sz,bcz), c1d.lapl2h(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, sz, bc, location = 'b')

def zblk(nx, nz, qx, qz, bc, location = 't'):
    """Create a block of zeros"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(nz,bcz),c1d.zblk(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, qx, qz, bc, location = location)

def i1j1(nx, nz, bc, coeff = 1.0):
    """Create operator for 1st integral in X and 1st integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i1(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 1, 1, bc)

def i1j1d1(nx, nz, bc, coeff = 1.0):
    """Create operator for 1st integral of 1st derivative in X and 1st integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.qid(nx,1,bcx))
    return c2dbc.constrain(mat, nx, nz, 1, 1, bc)

def i1j1e1(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 1st integral in X and 1st integral of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,1,bcz, zscale), c1d.i1(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 1, 1, bc)

def i2j1e1(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 2nd integral in X and 1st integral of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,1,bcz, zscale), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 1, bc)

def i2j2d1(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral of 1st derivative in X and 2nd integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i2d1(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2e1(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 2nd integral in X and 2nd integrazl of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(nz,bcz, zscale), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2d1e1(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 2nd integral of 1st derivative in X and 2nd integrazl of 1st derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(nz,bcz, zscale), c1d.i2d1(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2d2(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in X and 2nd integrazl of 2nd derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.qid(nx,2, bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2e2(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 2nd integral in X and 2nd integrazl of 2nd derivative in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz, zscale**2), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2d2e2(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create a quasi identity block of order 2,2"""

    return qid(nx,nz,2,2, bc, coeff*zscale**2)

def i2(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in X of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,0,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 0, bc)

def i2j1(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 1st integral in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 1, bc)

def i2j2(nx, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in X,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i2(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2laplh(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 2nd integral in x and 1st integral in z of horizontal Laplacian T_n(x)T_n(z)"""

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
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i2laplh(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 2, 2, bc)

def i2j2lapl(nx, nz, k, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 2nd integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(nz,bcz), c1d.i2laplh(nx,k,bcx))
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i2d2(nz,bcz,zscale**2), c1d.i2(nx,bcx))
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

def i4j1e1(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 4th integral in x and 1st integral of 1st derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,1,bcz,zscale), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j2e2(nx, nz, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 4th integral in x, 2nd integral of 2nd derivative in z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,2,bcz,zscale**2), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 2, bc)

def i4j4(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x,z of T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), c1d.i4(nx,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def i4j4d1(nx, nz, bc, coeff = 1.0):
    """Create operator for 4th integral of 1st derivative in X and 4th integral in Z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(nz,bcz), c1d.i4d1(nx,bcx))
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

def i4j4lapl(nx, nz, k, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 4th integral in x,z of Laplacian T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4laplh(nx,k,bcx))
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i4d2(nz,bcz,zscale**2), c1d.i4(nx,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def i4j1lapl2h(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 1st integral of z of horizontal Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i1(nz,bcz), c1d.i4lapl2h(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 1, bc)

def i4j2lapl2h(nx, nz, k, bc, coeff = 1.0):
    """Create operator for 4th integral in x and 2nd integral in z of horizontal Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(nz,bcz), c1d.i4lapl2h(nx,k,bcx))
    return c2dbc.constrain(mat, nx, nz, 4, 2, bc)

def i4j4lapl2(nx, nz, k, bc, coeff = 1.0, zscale = 1.0):
    """Create operator for 4th integral in x,z of Laplacian^2 T_n(x)T_n(z)"""

    bcx, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(nz,bcz), c1d.i4lapl2h(nx,k,bcx))
    bcx[0] = min(bcx[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + 2.0*spsp.kron(c1d.i4d2(nz,bcz,zscale**2), c1d.i4laplh(nx,k,bcx))
    mat = mat + spsp.kron(c1d.i4d4(nz,bcz,zscale**4), c1d.i4(nx,bcx))
    mat = coeff*mat
    return c2dbc.constrain(mat, nx, nz, 4, 4, bc)

def qid(nx, nz, qx, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qx in x in qz in z"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(nz,qz,bcz), c1d.qid(nx,qx,bcx))
    return c2dbc.constrain(mat, nx, nz, qx, qz, bc)

def sid(nx, nz, sx, sz, bc, coeff = 1.0):
    """Create a identity block order with last sx, sz rows zeroed"""

    bcx, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(nz,sz,bcz), c1d.sid(nx,sx,bcx))
    return c2dbc.constrain(mat, nx, nz, sx, sz, bc)
