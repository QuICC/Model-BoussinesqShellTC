"""Module provides functions to generate sparse operators in a cartesian box."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.geometry.cartesian.cartesian_boundary_3d as c3dbc


def convert_bc(bc):
    """Convert boundary dictionary into x, y and z kronecker product boundaries"""

    if bc['x'][0] < 0:
        bcx = bc['x']
    else:
        bcx = c1d.c1dbc.no_bc()
        for key, val in bc['x'].items():
            if key != 0:
                bcx[key] = val

    if bc['y'][0] < 0:
        bcy = bc['y']
    else:
        bcy = c1d.c1dbc.no_bc()
        for key, val in bc['y'].items():
            if key != 0:
                bcy[key] = val

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = c1d.c1dbc.no_bc()
        for key, val in bc['z'].items():
            if key != 0:
                bcz[key] = val

    return (bcx, bcy, bcz)

def d1(nx, ny, nz, bc, coeff = 1.0, sx = 1, sy = 0, sz = 0):
    """Create operator for 1st X derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(ny, sy, bcy), c2d.d1(nx,nz,{'x':bcx,'z':bcz}, sx = sx, sz = sz))
    return  c3dbc.constrain(mat, nx, ny, nz, sx, sy, sz, bc, location = 'b')

def e1(nx, ny, nz, bc, coeff = 1.0, sx = 0, sy = 1, sz = 0):
    """Create operator for 1st Y derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.d1(ny,bcy, zr = sy), c2d.d1(nx,nz,{'x':bcx,'z':bcz}, sx = sx, sz = sz))
    return  c3dbc.constrain(mat, nx, ny, nz, sx, sy, sz, bc, location = 'b')

def f1(nx, ny, nz, bc, coeff = 1.0, sx = 0, sy = 0, sz = 1):
    """Create operator for 1st Z derivative of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.sid(ny, sy, bcy), c2d.e1(nx,nz,{'x':bcx,'z':bcz}, sx = sx, sz = sz))
    return  c3dbc.constrain(mat, nx, ny, nz, sx, sy, sz, bc, location = 'b')

def zblk(nx, ny, nz, qx, qy, qz, bc, location = 't'):
    """Create a block of zeros"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.zblk(ny,bcy),spsp.kron(c1d.zblk(nz,bcz),c1d.zblk(nx,bcx)))
    return c3dbc.constrain(mat, nx, ny, nz, qx, qy, qz, bc, location = location)

def i2j2k2d2d2d2(nx, ny, nz, bc, coeff = 1.0):
    """Create a quasi identity block of order 2,2,2"""

    return qid(nx,ny,nz,2,2,2,bc,coeff)

def i2j2k2(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(ny,bcy), c2d.i2j2(nx,nz,{'x':bcx,'z':bcz}))
    return  c3dbc.constrain(mat, nx, ny, nz, 2, 2, 2, bc)

def i2j2k2d1(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,y,z of Laplacian T'_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(ny,bcy), c2d.i2j2d1(nx,nz,{'x':bcx,'z':bcz}))
    return  c3dbc.constrain(mat, nx, ny, nz, 2, 2, 2, bc)

def i2j2k2e1(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T'_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2d1(ny,bcy), c2d.i2j2(nx,nz,{'x':bcx,'z':bcz}))
    return  c3dbc.constrain(mat, nx, ny, nz, 2, 2, 2, bc)

def i2j2k2f1(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T'_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i2(ny,bcy), c2d.i2j2e1(nx,nz,{'x':bcx,'z':bcz}))
    return  c3dbc.constrain(mat, nx, ny, nz, 2, 2, 2, bc)

def i2j2k2lapl(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i2(ny,bcy),spsp.kron(c1d.i2(nz,bcz),c1d.i2d2(nx,bcx)))
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i2d2(ny,bcy),spsp.kron(c1d.i2(nz,bcz),c1d.i2(nx,bcx)))
    mat = mat + spsp.kron(c1d.i2(ny,bcy),spsp.kron(c1d.i2d2(nz,bcz),c1d.i2(nx,bcx)))
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, ny, nz, 2, 2, 2, bc)

def i4j4k4(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.i4(ny,bcy), c2d.i4j4(nx,nz,{'x':bcx,'z':bcz}))
    return  c3dbc.constrain(mat, nx, ny, nz, 4, 4, 4, bc)

def i4j4k4lapl(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d2(nx,bcx)))
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4(nx,bcx)))
    mat = mat + spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4(nx,bcx)))
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, ny, nz, 4, 4, 4, bc)

def i4j4k4lapl2(nx, ny, nz, bc, coeff = 1.0):
    """Create operator for 4th integral in x,y,z of Laplacian^2 T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d4(nx,bcx)))
    bcx[0] = min(bcx[0], 0)
    bcy[0] = min(bcy[0], 0)
    bcz[0] = min(bcz[0], 0)
    mat = mat + spsp.kron(c1d.i4d4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4(nx,bcx)))
    mat = mat + spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d4(nz,bcz),c1d.i4(nx,bcx)))
    mat = mat + 2*spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d2(nx,bcx)))
    mat = mat + 2*spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4(nx,bcx)))
    mat = mat + 2*spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4d2(nx,bcx)))
    mat = coeff*mat
    return  c3dbc.constrain(mat, nx, ny, nz, 4, 4, 4, bc)

def qid(nx, ny, nz, qx, qy, qz, bc, coeff = 1.0):
    """Create a quasi identity block order qx,qy,qz"""

    bcx, bcy, bcz = convert_bc(bc)
    mat = coeff*spsp.kron(c1d.qid(ny,qy,bcy), c2d.qid(nx,nz,qx,qz,{'x':bcx,'z':bcz}))
    return c3dbc.constrain(mat, nx, ny, nz, qx, qy, qz, bc)
