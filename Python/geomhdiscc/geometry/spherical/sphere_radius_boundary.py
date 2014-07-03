"""Module provides functions to generate the radial boundary conditions in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np


def constrain(mat, l, bc, eq_zrows):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, l, bc)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, l, bc, eq_zrows)
    else:
        bc_mat = mat

    return bc_mat


def apply_tau(mat, l, bc):
    """Add Tau lines to the matrix"""

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], l%2, bc[1:])
    if bc[0] == 11:
        cond = tau_diff(mat.shape[0], l%2, bc[1:])
    elif bc[0] == 20:
        cond = tau_value_diff(mat.shape[0], l%2, bc[1:])
    elif bc[0] == 21:
        cond = tau_value_diff2(mat.shape[0], l%2, bc[1:])

    if cond.dtype == 'complex_':
        bc_mat = mat.astype('complex_').tolil()
    else:
        bc_mat = mat.tolil()

    bc_mat[0:cond.shape[0],:] = cond

    return bc_mat


def tau_value(nr, parity, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    cond.append([c*norm_c(i) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)


def tau_diff(nr, parity, coeffs = None):
    """Create the first derivative tau line(s)"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    cond.append([c*i**2 for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)


def tau_diff2(nr, parity, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    cond.append([c*((1/3)*(i**4 - i**2)) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)


def tau_value_diff(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    cond.append(list(tau_value(nr,parity,coeffs)[0]))
    cond.append(list(tau_diff(nr,parity,coeffs)[0]))

    return np.array(cond)


def tau_value_diff2(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    cond.append(list(tau_value(nr,parity,coeffs)[0]))
    cond.append(list(tau_diff2(nr,parity,coeffs)[0]))

    return np.array(cond)


def apply_galerkin(mat, l, bc, eq_zero_rows):
    """Apply a Galerkin stencil on the matrix"""

    return mat


def norm_c(n):
    """Compute the chebyshev normalisation c factor"""

    if n > 0:
        return 2
    else:
        return 1
