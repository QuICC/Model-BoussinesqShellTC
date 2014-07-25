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

def qid(nx, q, bc):
    """Create a quasi indentity"""

    if bc[0] < 0:
        mat = spsp.identity(nx-bc[0]//10)
    else:
        offsets = [0]
        diags = [[0]*q + [1]*(nx-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def bid(nx, q, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.identity(nx-bc[0]//10)
    else:
        offsets = [-q]
        diags = [[1]*(nx-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def constrain(mat, nr, nz, qr, qz, bc):
    """Contrain the matrix with the Tau boundary condition"""

    sr = 0
    sz = qz

    bc_mat = mat
    if bc['r'][0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, bc['r'])
        bc_mat = bc_mat + spsp.kron(qid(nz,sz,bc['z']), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'])
        bc_mat = bc_mat + spsp.kron(bcMat, qid(nr,sr,bc['r']))

    return bc_mat
