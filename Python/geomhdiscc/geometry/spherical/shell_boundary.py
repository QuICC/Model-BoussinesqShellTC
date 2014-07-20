"""Module provides functions to generate the boundary conditions in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.spherical.shell_radius_boundary as radbc
from geomhdiscc.geometry.spherical.shell_radius_boundary import no_bc


def qid(n, q, bc):
    """Create a quasi indentity"""

    if bc[0] < 0:
        mat = spsp.identity(n-bc[0]//10)
    else:
        offsets = [0]
        diags = [[0]*q + [1]*(n-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def bid(n, q, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.identity(n-bc[0]//10)
    else:
        offsets = [-q]
        diags = [[1]*(n-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()

def constrain(mat, nr, maxl, m, bc):
    """Contrain the matrix with the tau boundary condition"""

    bc_mat = mat
    if bc[0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bc_mat = radbc.constrain(bcMat, bc)
        for l in range(m+1, maxl+1):
            bcMat = spsp.lil_matrix((nr,nr))
            bcMat = radbc.constrain(bcMat, bc)
            bc_mat = spsp.block_diag((bc_mat,bcMat))

        bc_mat = mat + bc_mat

    return bc_mat
