"""Module provides functions to generate the radial boundary conditions in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import itertools

import geomhdiscc.base.utils as utils


use_parity_bc = False

def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, bc, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, bc, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, bc)
    else:
        bc_mat = mat

    # top row(s) restriction if required
    if bc.get('rt', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], 'rt', bc['rt'])*bc_mat

    # bottom row(s) restriction if required
    if bc.get('rb', 0) > 0:
        bc_mat = restrict_eye(bc_mat.shape[0], 'rb', bc['rb'])*bc_mat

    # left columns restriction if required
    if bc.get('cl', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], 'cl', bc['cl'])

    # right columns restriction if required
    if bc.get('cr', 0) > 0:
        bc_mat = bc_mat*restrict_eye(bc_mat.shape[1], 'cr', bc['cr'])

    # top row(s) zeroing if required
    if bc.get('zt', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[0:bc['zt'],:] = 0
        bc_mat = bc_mat.tocoo()

    # bottom row(s) zeroing if required
    if bc.get('zb', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[-bc['zb']:,:] = 0
        bc_mat = bc_mat.tocoo()

    # left columns zeroing if required
    if bc.get('zl', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[:, 0:bc['zt']] = 0
        bc_mat = bc_mat.tocoo()

    # right columns zeroing if required
    if bc.get('zr', 0) > 0:
        bc_mat = bc_mat.tolil()
        bc_mat[:, -bc['zr']:] = 0
        bc_mat = bc_mat.tocoo()

    return bc_mat

def apply_tau(mat, bc, location = 't'):
    """Add Tau lines to the matrix"""

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], 1, bc.get('c',None))
    elif bc[0] == 11:
        cond = tau_value(mat.shape[0], -1, bc.get('c',None))
    elif bc[0] == 12:
        cond = tau_diff(mat.shape[0], 1, bc.get('c',None))
    elif bc[0] == 13:
        cond = tau_diff(mat.shape[0], -1, bc.get('c',None))
    elif bc[0] == 20:
        cond = tau_value(mat.shape[0], 0, bc.get('c',None))
    elif bc[0] == 21:
        cond = tau_diff(mat.shape[0], 0, bc.get('c',None))
    elif bc[0] == 40:
        cond = tau_value_diff(mat.shape[0], 0, bc.get('c',None))
    elif bc[0] == 41:
        cond = tau_value_diff2(mat.shape[0], 0, bc.get('c',None))

    if cond.dtype == 'complex_':
        bc_mat = mat.astype('complex_').tolil()
    else:
        bc_mat = mat.tolil()

    if location == 't':
        bc_mat[0:cond.shape[0],:] = cond
    elif location == 'b':
        bc_mat[-cond.shape[0]:,:] = cond

    return bc_mat

def tau_value(nr, pos, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None:
        it = itertools.cycle([1.0])
    else:
        try:
            if len(coeffs) == (1 + (pos == 0)):
                it = iter(coeffs)
            elif len(coeffs) == 1:
                it = itertools.cycle(coeffs)
            else:
                raise RuntimeError
        except:
            it = itertools.cycle([coeffs])

    cond = []
    c = next(it)
    if pos >= 0:
        cond.append([c*tau_c(i) for i in np.arange(0,nr)])
        c = next(it)

    if pos <= 0:
        cond.append([c*tau_c(i)*(-1.0)**i for i in np.arange(0,nr)])

    if use_parity_bc and pos == 0:
        t = cond[0]
        cond[0] = [(cond[0][i] + cond[1][i])/2 for i in np.arange(0,nx)]
        cond[1] = [(t[i] - cond[1][i])/2 for i in np.arange(0,nx)]

    return np.array(cond)

def tau_diff(nr, pos, coeffs = None):
    """Create the first derivative tau line(s)"""

    if coeffs is None:
        it = itertools.cycle([1.0])
    else:
        try:
            if len(coeffs) == (1 + (pos == 0)):
                it = iter(coeffs)
            elif len(coeffs) == 1:
                it = itertools.cycle(coeffs)
            else:
                raise RuntimeError
        except:
            it = itertools.cycle([coeffs])

    cond = []
    c = next(it)
    if pos >= 0:
        cond.append([c*i**2 for i in np.arange(0,nr)])
        c = next(it)

    if pos <= 0:
        cond.append([(-1.0)**(i+1)*c*i**2 for i in np.arange(0,nr)])

    if use_parity_bc and pos == 0:
        t = cond[0]
        cond[0] = [(cond[0][i] + cond[1][i])/2 for i in np.arange(0,nx)]
        cond[1] = [(t[i] - cond[1][i])/2 for i in np.arange(0,nx)]

    return np.array(cond)

def tau_diff2(nr, pos, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        it = itertools.cycle([1.0])
    else:
        try:
            if len(coeffs) == (1 + (pos == 0)):
                it = iter(coeffs)
            elif len(coeffs) == 1:
                it = itertools.cycle(coeffs)
            else:
                raise RuntimeError
        except:
            it = itertools.cycle([coeffs])

    cond = []
    c = next(it)
    if pos >= 0:
        cond.append([c*((1/3)*(i**4 - i**2)) for i in np.arange(0,nr)])
        c = next(it)

    if pos <= 0:
        cond.append([c*(((-1.0)**i/3)*(i**4 - i**2)) for i in np.arange(0,nr)])

    if use_parity_bc and pos == 0:
        t = cond[0]
        cond[0] = [(cond[0][i] + cond[1][i])/2 for i in np.arange(0,nx)]
        cond[1] = [(t[i] - cond[1][i])/2 for i in np.arange(0,nx)]

    return np.array(cond)

def tau_value_diff(nr, pos, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    if pos >= 0:
        cond.append(list(tau_value(nr,1,coeffs)[0]))
        cond.append(list(tau_diff(nr,1,coeffs)[0]))

    if pos <= 0:
        cond.append(list(tau_value(nr,-1,coeffs)[0]))
        cond.append(list(tau_diff(nr,-1,coeffs)[0]))

    if use_parity_bc and pos == 0:
        tv = cond[0]
        td = cond[1]
        cond[0] = [(cond[0][i] + cond[2][i])/2 for i in np.arange(0,nx)]
        cond[1] = [(cond[1][i] + cond[3][i])/2 for i in np.arange(0,nx)]
        cond[2] = [(tv[i] - cond[2][i])/2 for i in np.arange(0,nx)]
        cond[3] = [(td[i] - cond[3][i])/2 for i in np.arange(0,nx)]

    return np.array(cond)

def tau_value_diff2(nr, pos, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    if pos >= 0:
        cond.append(list(tau_value(nr,1,coeffs)[0]))
        cond.append(list(tau_diff2(nr,1,coeffs)[0]))

    if pos <= 0:
        cond.append(list(tau_value(nr,-1,coeffs)[0]))
        cond.append(list(tau_diff2(nr,-1,coeffs)[0]))

    if use_parity_bc and pos == 0:
        tv = cond[0]
        td = cond[1]
        cond[0] = [(cond[0][i] + cond[2][i])/2 for i in np.arange(0,nx)]
        cond[1] = [(cond[1][i] + cond[3][i])/2 for i in np.arange(0,nx)]
        cond[2] = [(tv[i] - cond[2][i])/2 for i in np.arange(0,nx)]
        cond[3] = [(td[i] - cond[3][i])/2 for i in np.arange(0,nx)]

    return np.array(cond)

def apply_galerkin(mat, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = mat*stencil(nr, parity, bc)
    return mat

def restrict_eye(nr, t, q):
    """Create the non-square identity to restrict matrix"""

    if t == 'rt':
        offsets = [q]
        diags = [[1]*(nr-q)]
        nrows = nr - q
        ncols = nr
    elif t == 'rb':
        offsets = [0]
        diags = [[1]*(nr-q)]
        nrows = nr - q
        ncols = nr
    elif t == 'cl':
        offsets = [-q]
        diags = [[1]*(nr-q)]
        nrows = nr
        ncols = nr - q
    elif t == 'cr':
        offsets = [0]
        diags = [[1]*(nr-q)]
        nrows = nr
        ncols = nr - q

    return spsp.diags(diags, offsets, (nrows, ncols))

def tau_c(n):
    """Compute the chebyshev normalisation c factor"""
    if n > 0:
        return 2
    else:
        return 1

def galerkin_c(n):
    """Compute the chebyshev normalisation c factor for galerkin boundary"""

    if n > 0:
        return 1
    else:
        return 0.5
