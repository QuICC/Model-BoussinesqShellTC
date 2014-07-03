"""Module provides functions to generate the boundary conditions in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.geometry.spherical.shell_radius_boundary as radbc

def qid(nr, q, bc):
    """Create a quasi indentity"""

    if bc[0] < 0:
        mat = spsp.identity(nr-bc[0]//10)
    else:
        offsets = [0]
        diags = [[0]*q + [1]*(nr-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()


def bid(nr, q, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.identity(nr-bc[0]//10)
    else:
        offsets = [-q]
        diags = [[1]*(nr-q)]

        mat = spsp.diags(diags, offsets)

    return mat.tocsr()


def constrain(mat, nr, nl, bc):
    """Contrain the matrix with the tau boundary condition"""

    bc_mat = mat
    if bc[0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, bc, 0)
        bc_mat = bc_mat + spsp.kron(bid(nl,0,[0]), bcMat)

    return bc_mat
