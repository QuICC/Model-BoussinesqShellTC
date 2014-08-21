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
    tbc['cr'] = bc.get('cr', 0)
    tbc['rb'] = bc.get('rb', 0)
    tbc['cl'] = bc.get('cl', 0)
    tbc['rt'] = bc.get('rt', 0)
    return c1dbc.constrain(mat, tbc)

def constrain(mat, nx, ny, nz, qx, qy, qz, bc, location = 't'):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'x')
    if priority == 'x':
        sx = qx
        dx = 0
        sy = 0
        dy = 0
        sz = 0
        dz = 0
    elif priority == 'y':
        sx = 0
        dx = 0
        sy = qy
        dy = 0
        sz = 0
        dz = 0
    elif priority == 'z':
        sx = 0
        dx = 0
        sy = 0
        dy = 0
        sz = qz
        dz = 0
    elif priority == 'n':
        sx = qx
        dx = 0
        sy = qy
        dy = 0
        sz = qz
        dz = 0
    if priority == 'sx':
        sx = qx
        dx = -qx
        sy = 0
        dy = 0
        sz = 0
        dz = 0
    elif priority == 'sy':
        sx = 0
        dx = 0
        sy = qy
        dy = -qy
        sz = 0
        dz = 0
    elif priority == 'sz':
        sx = 0
        dx = 0
        sy = 0
        dy = 0
        sz = qz
        dz = -qz
            
    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, bc['x'], location = location)
        bc_mat = bc_mat + spsp.kron(bid(nz,sz,dz,bc['z'], location = location), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], location = location)
        if bc['x'][0] >= 0:
            bc_mat = bc_mat + spsp.kron(bcMat, bid(nx,sx,dx,bc['x'], location = location))

        else:
            tmpB = c1dbc.constrain(bid(nx,0,0,c1dbc.no_bc(), location = location),bc['x'], location = location)
            bc_mat = bc_mat + spsp.kron(bcMat, tmpB)

    return bc_mat
