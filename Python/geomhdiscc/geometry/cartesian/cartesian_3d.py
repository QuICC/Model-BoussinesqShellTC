"""Module provides functions to generate sparse operators in a cartesian box."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d

def convert_bc(bc):
    """Convert boundary dictionary into x, y and z kronecker product boundaries"""

    if bc['x'][0] < 0:
        bcx = bc['x']
    else:
        bcx = [0]

    if bc['y'][0] < 0:
        bcy = bc['y']
    else:
        bcy = [0]

    if bc['z'][0] < 0:
        bcz = bc['z']
    else:
        bcz = [0]

    return (bcx, bcy, bcz)


def zblk(nx, ny, nz, qx, qy, qz, bc):
    """Create a block of zeros"""

    bcx, bcy, bcz = convert_bc(bc)
    return spsp.kron(c1d.zblk(ny,qy,bcy),spsp.kron(c1d.zblk(nz,qz,bcz),c1d.zblk(nx,qx,bcx)))


def i2j2k2d2d2d2(nx, ny, nz, bc):
    """Create a quasi identity block of order 2,2,2"""

    return qid(nx,ny,nz,2,2,2,bc)


def i2j2k2(nx, ny, nz, bc):
    """Create operator for 2nd integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    return spsp.kron(c1d.i2(ny,bcy), c2d.i2j2(nx,nz,{'x':bcx,'z':bcz}))


def i2j2k2lapl(nx, ny, nz, bc):
    """Create operator for 2nd integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    op = spsp.kron(c1d.i2(ny,bcy),spsp.kron(c1d.i2(nz,bcz),c1d.i2d2(nx,bcx)))
    op = op + spsp.kron(c1d.i2d2(ny,bcy),spsp.kron(c1d.i2(nz,bcz),c1d.i2(nx,bcx)))
    op = op + spsp.kron(c1d.i2(ny,bcy),spsp.kron(c1d.i2d2(nz,bcz),c1d.i2(nx,bcx)))

    return  op


def i4j4k4(nx, ny, nz, bc):
    """Create operator for 4th integral in x,y,z of T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    return spsp.kron(c1d.i4(ny,bcy), c2d.i4j4(nx,nz,{'x':bcx,'z':bcz}))


def i4j4k4lapl(nx, ny, nz, bc):
    """Create operator for 4th integral in x,y,z of Laplacian T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    op = spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d2(nx,bcx)))
    op = op + spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4(nx,bcx)))
    op = op + spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4(nx,bcx)))

    return  op


def i4j4k4lapl2(nx, ny, nz, bc):
    """Create operator for 4th integral in x,y,z of Laplacian^2 T_n(x)T_n(y)T_n(z)"""

    bcx, bcy, bcz = convert_bc(bc)
    op = spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d4(nx,bcx)))
    op = op + spsp.kron(c1d.i4d4(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4(nx,bcx)))
    op = op + spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d4(nz,bcz),c1d.i4(nx,bcx)))
    op = op + 2*spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4(nz,bcz),c1d.i4d2(nx,bcx)))
    op = op + 2*spsp.kron(c1d.i4d2(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4(nx,bcx)))
    op = op + 2*spsp.kron(c1d.i4(ny,bcy),spsp.kron(c1d.i4d2(nz,bcz),c1d.i4d2(nx,bcx)))

    return  op


def qid(nx, ny, nz, qx, qy, qz, bc):
    """Create a quasi identity block order qx,qy,qz"""

    bcx, bcy, bcz = convert_bc(bc)
    return spsp.kron(c1d.qid(ny,qy,bcy), c2d.qid(nx,nz,qx,qz,{'x':bcx,'z':bcz}))
