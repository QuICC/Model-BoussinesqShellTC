"""Module provides functions to generate the radial boundary conditions in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, l, bc, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, l, bc, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, l, bc)
    else:
        bc_mat = mat

    # Restrict if required
    if bc.get('r', 0) > 0:
        assert bc[0] <= 0
        bc_mat = stencil_eye(mat.shape[0], bc['r'])*bc_mat

    return bc_mat

def apply_tau(mat, l, bc, location = 't'):
    """Add Tau lines to the matrix"""

    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], l%2, bc.get('c',None))
    if bc[0] == 11:
        cond = tau_diff(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 20:
        cond = tau_value_diff(mat.shape[0], l%2, bc.get('c',None))
    elif bc[0] == 21:
        cond = tau_value_diff2(mat.shape[0], l%2, bc.get('c',None))
    # Set last modes to zero
    elif bc[0] > 990 and bc[0] < 1000:
        cond = tau_last(mat.shape[1], bc[0]-990)
        nbc = bc[0]-990

    if not spsp.isspmatrix_coo(mat):
        mat = mat.tocoo()
    if location == 't':
        s = 0
    elif location == 'b':
        s = mat.shape[0]-nbc

    conc = np.concatenate
    for i,c in enumerate(cond):
        mat.data = conc((mat.data, c))
        mat.row = conc((mat.row, [s+i]*mat.shape[1]))
        mat.col = conc((mat.col, np.arange(0,mat.shape[1])))

    return mat

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
    cond.append([c*((1.0/3.0)*(i**4 - i**2)) for i in np.arange(parity, 2*nr, 2)])

    return np.array(cond)

def tau_value_diff(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    cond.append(list(tau_value(nr,parity,coeffs)[0]))
    cond.append(list(tau_diff(nr,parity,coeffs)[0]))

    return np.array(cond)

def tau_value_diff2(nr, parity, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    cond.append(list(tau_value(nr,parity,coeffs)[0]))
    cond.append(list(tau_diff2(nr,parity,coeffs)[0]))

    return np.array(cond)

def apply_galerkin(mat, l, bc):
    """Apply a Galerkin stencil on the matrix"""

    return mat

def norm_c(n):
    """Compute the chebyshev normalisation c factor"""

    if n > 0:
        return 2
    else:
        return 1
