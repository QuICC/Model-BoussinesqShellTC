"""Module provides functions to generate the boundary conditions in a sphere with Chebyshev expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.spherical.sphere_radius_boundary_chebyshev as radbc
from geomhdiscc.geometry.spherical.sphere_radius_boundary_chebyshev import no_bc


def no_bc():
    """Get a no boundary condition flag"""

    return radbc.no_bc()

def ldependent_bc(bc, l):
    """Set harmonic degree l for boundary conditions"""

    if bc.get('c', None) is not None:
        if bc['c'].get('l', None) is not None:
            bc['c']['l'] = l

    return bc

def constrain(mat, nr, maxnl, m, bc, zero_l_zero = False, restriction = None):
    """Contrain the matrix with the tau boundary condition"""

    if bc[0] > 0:
        bcMat = spsp.coo_matrix((nr,nr))
        if m == 0 and zero_l_zero:
            bc_mat = bcMat
        else:
            if restriction is None or m in restriction:
                bc = ldependent_bc(bc, m)
                bc_mat = radbc.constrain(bcMat, m, bc)
            else:
                bc_mat = bcMat
        for l in range(m+1, maxnl):
            bcMat = spsp.coo_matrix((nr,nr))
            if restriction is None or l in restriction:
                bc = ldependent_bc(bc, l)
                bcMat = radbc.constrain(bcMat, l, bc)
            bc_mat = spsp.block_diag((bc_mat,bcMat), format = 'coo')

        bc_mat = mat + bc_mat

        return bc_mat
    else:
        return mat