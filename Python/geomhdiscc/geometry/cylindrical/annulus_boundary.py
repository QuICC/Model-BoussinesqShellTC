"""Module provides functions to generate the boundary conditions in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import geomhdiscc.geometry.cylindrical.annulus_radius_boundary as radbc

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


def constrain(mat, nr, nz, bc, eq_zrows_r, eq_zrows_z):
    """Contrain the matrix with the Tau boundary condition"""

    bc_mat = mat
    if bc['r'][0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, bc['r'], 0)
        bc_mat = bc_mat + spsp.kron(bid(nz,eq_zrows_z,bc['z']), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], 0)
        bc_mat = bc_mat + spsp.kron(bcMat, bid(nr,0,bc['x']))

    return bc_mat
