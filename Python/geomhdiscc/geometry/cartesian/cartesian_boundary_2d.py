"""Module provides functions to generate the boundary conditions in a cartesian 2D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc


def no_bc():
    """Get a no boundary condition flag"""

    return {'x':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

def bid(nx, q, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.identity(nx-bc['r'])
    else:
        offsets = [-q]
        diags = [[1]*(nx-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def qid(n, q, bc):
    """Create a quasi indentity"""

    if bc[0] < 0:
        mat = spsp.identity(n-bc[0]//10)
    else:
        offsets = [0]
        diags = [[0]*q + [1]*(n-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def constrain(mat, nx, nz, bc):
    """Contrain the matrix with the Tau boundary condition"""

    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, bc['x'])
        bc_mat = bc_mat + spsp.kron(bid(nz,bc['z'][0]//10,bc['z']), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'])
        bc_mat = bc_mat + spsp.kron(bcMat, bid(nx,0*bc['x'][0]//10,bc['x']))

    return bc_mat
