"""Module provides functions to generate the boundary conditions in a cartesian 3D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.polynomial.chebyshev as cheby
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import geomhdiscc.transform.cartesian as phys


def no_bc():
    """Get a no boundary condition flag"""

    return {'x':c1dbc.no_bc(), 'y':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

def bid(nx, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(nx-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = c1dbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return c1dbc.constrain(mat, tbc)

def constrain(mat, nx, ny, nz, qx, qy, qz, bc, location = 't'):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'xz')
    sx = 2*[0]
    dx = 2*[0]
    sy = 2*[0]
    dy = 2*[0]
    sz = 2*[0]
    dz = 2*[0]
    if priority == 'xy':
        sx = [qx,qx]
        sy = [0,qy]
    elif priority == 'xz':
        sx = [qx,qx]
        sz = [0,qz]
    elif priority == 'yx':
        sx = [0,qx]
        sy = [qy,qy]
    elif priority == 'yz':
        sy = [qy,qy]
        sz = [qz,0]
    elif priority == 'zx':
        sx = [qx,0]
        sz = [qz,qz]
    elif priority == 'zy':
        sy = [qy,0]
        sz = [qz,qz]
    elif priority == 'n':
        sx = [qx,qx]
        sy = [qy,qy]
        sz = [qz,qz]
            
    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, bc['x'], location = location)
        bc_mat = bc_mat + spsp.kron(bid(ny,sy[0],dy[0],bc['y'], location = location), spsp.kron(bid(nz,sz[0],dz[0],bc['z']),bcMat))

    if bc['y'][0] > 0:
        bcMat = spsp.lil_matrix((ny,ny))
        bcMat = c1dbc.constrain(bcMat, bc['y'], location = location)
        bc_mat = bc_mat + spsp.kron(bcMat, spsp.kron(bid(nz,sz[1],dz[1],bc['z']),bid(nx,sx[0],dx[0],bc['x'], location = location)))

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], location = location)
        bc_mat = bc_mat + spsp.kron(bid(ny,sy[1],dy[1],bc['y'], location = location), spsp.kron(bcMat,bid(nx,sx[1],dx[1],bc['x'])))

    return bc_mat