"""Module provides functions to generate the boundary conditions in a cartesian 2D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc


def no_bc():
    """Get a no boundary condition flag"""

    return {'x':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

def qid(n, q, bc):
    """Create a quasi indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [0]
        diags = [[0]*q + [1]*(n-q)]

        mat = spsp.diags(diags, offsets)

    tbc = c1dbc.no_bc()
    tbc['cr'] = bc.get('cr', 0)
    tbc['rb'] = bc.get('rb', 0)
    tbc['cl'] = bc.get('cl', 0)
    tbc['rt'] = bc.get('rt', 0)
    return c1dbc.constrain(mat, tbc)

def bid(nx, q, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [-q]
        diags = [[1]*(nx-q)]

        mat = spsp.diags(diags, offsets)

    tbc = c1dbc.no_bc()
    tbc['cr'] = bc.get('cr', 0)
    tbc['rb'] = bc.get('rb', 0)
    tbc['cl'] = bc.get('cl', 0)
    tbc['rt'] = bc.get('rt', 0)
    return c1dbc.constrain(mat, tbc)

def constrain(mat, nx, nz, qx, qz, bc):
    """Contrain the matrix with the Tau boundary condition"""

    sx = 0
    sz = qz
            
    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, bc['x'])
        bc_mat = bc_mat + spsp.kron(qid(nz,sz,bc['z']), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'])
        if bc['x'][0] >= 0:
            bc_mat = bc_mat + spsp.kron(bcMat, qid(nx,sx,bc['x']))
            #bc_mat = bc_mat + spsp.kron(bcMat, bid(nx,sx,bc['x']))
        else:
            tmpB = c1dbc.constrain(qid(nx,0,c1dbc.no_bc()),bc['x'])
            bc_mat = bc_mat + spsp.kron(bcMat, tmpB)

    return bc_mat
