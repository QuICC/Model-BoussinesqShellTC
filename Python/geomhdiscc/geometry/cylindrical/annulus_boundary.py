"""Module provides functions to generate the boundary conditions in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import geomhdiscc.geometry.cylindrical.annulus_radius_boundary as radbc


def no_bc():
    """Get a no boundary condition flag"""

    return {'r':radbc.no_bc(), 'z':c1dbc.no_bc()}

def bid(n, q, d, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [-d]
        diags = [[1]*(n-d)]

        mat = spsp.diags(diags, offsets).tolil()
        mat[0:q,:] = 0

    return mat.tocsr()

def constrain(mat, nr, nz, qr, qz, bc, location = 't'):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'r')
    if priority == 'r':
        sr = qr
        sz = 0
    elif priority == 'z':
        sr = 0
        sz = qz
    elif priority == 'n':
        sr = qr
        sz = qz

    bc_mat = mat
    if bc['r'][0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, bc['r'], location = location)
        bc_mat = bc_mat + spsp.kron(bid(nz, sz, 0, bc['z'], location = location), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], location = location)
        bc_mat = bc_mat + spsp.kron(bcMat, bid(nr, sr, 0, bc['r'], location = location))

    return bc_mat
