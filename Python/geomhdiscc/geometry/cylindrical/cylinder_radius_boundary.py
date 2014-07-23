"""Module provides functions to generate the radial boundary conditions in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, m, bc):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, m, bc)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, m, bc)
    else:
        bc_mat = mat

    # Restrict if required
    if bc.get('r', 0) > 0:
        assert bc[0] <= 0
        bc_mat = stencil_eye(mat.shape[0], bc['r'])*bc_mat

    return bc_mat

def apply_tau(mat, m, bc):
    """Add Tau lines to the matrix"""

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], m%2, bc.get('c',None))
    elif bc[0] == 11:
        cond = tau_diff(mat.shape[0], m%2, bc.get('c',None))
    elif bc[0] == 13:
        cond = tau_1rdr(mat.shape[0], m%2, bc.get('c',None))

    if cond.dtype == 'complex_':
        bc_mat = mat.astype('complex_').tolil()
    else:
        bc_mat = mat.tolil()

    bc_mat[0:cond.shape[0],:] = cond

    return bc_mat

def tau_value(nr, parity, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*norm_c(i) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_diff(nr, parity, coeffs = None):
    """Create the first derivative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*i**2 for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_diff2(nr, parity, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append([c*((1/3)*(i**4 - i**2)) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_1rdr(nr, parity, coeffs = None):
    """Create the 1/r D_r tau line(s)"""

    cond = tau_value(nr, parity, coeffs) + tau_diff(nr, parity, coeffs)

    return cond

def apply_galerkin(mat, m, bc):
    """Apply a Galerkin stencil on the matrix"""

    return mat

def norm_c(n):
    """Compute the chebyshev normalisation c factor"""

    if n > 0:
        return 2
    else:
        return 1
