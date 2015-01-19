"""Module provides functions to generate the boundary conditions in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.spherical.shell_radius_boundary as radbc
from geomhdiscc.geometry.spherical.shell_radius_boundary import no_bc


def no_bc():
    """Get a no boundary condition flag"""

    return radbc.no_bc()

def constrain(mat, nr, maxnl, m, bc, zero_l_zero = False):
    """Contrain the matrix with the tau boundary condition"""

    bc_mat = mat
    if bc[0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        if m == 0 and zero_l_zero:
            bc_mat = bcMat
        else:
            bc_mat = radbc.constrain(bcMat, bc)
        for l in range(m+1, maxnl):
            bcMat = spsp.lil_matrix((nr,nr))
            bcMat = radbc.constrain(bcMat, bc)
            bc_mat = spsp.block_diag((bc_mat,bcMat))

        bc_mat = mat + bc_mat

    return bc_mat
