"""Module provides functions to generate the radial boundary conditions in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import itertools


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

    # Restrict if required
    if bc.get('r', 0) > 0:
        assert bc[0] <= 0
        bc_mat = stencil_eye(mat.shape[0], bc['r'])*bc_mat

    return bc_mat

def apply_tau(mat, bc, location = 't'):
    """Add Tau lines to the matrix"""

    # Outer u = 0
    if bc[0] == 10:
        cond = tau_value(mat.shape[0], 1, bc.get('c',None))
    # Inner u = 0
    elif bc[0] == 11:
        cond = tau_value(mat.shape[0], -1, bc.get('c',None))
    # Outer D u = 0
    elif bc[0] == 12:
        cond = tau_diff(mat.shape[0], 1, bc.get('c',None))
    # Inner D u = 0
    elif bc[0] == 13:
        cond = tau_diff(mat.shape[0], -1, bc.get('c',None))
    # integral
    elif bc[0] == 14:
        cond = tau_integral(mat.shape[0], 1, bc.get('c',None))
    # u = 0
    elif bc[0] == 20:
        cond = tau_value(mat.shape[0], 0, bc.get('c',None))
    # D u = 0
    elif bc[0] == 21:
        cond = tau_diff(mat.shape[0], 0, bc.get('c',None))
    # D^2 u = 0
    elif bc[0] == 22:
        cond = tau_diff2(mat.shape[0], 0, bc.get('c',None))
    # 1/r D(r u) = 0
    elif bc[0] == 23:
        cond = tau_1rdr(mat.shape[0], 0, bc.get('c',None))
    # 1/r D(u/r) = 0
    elif bc[0] == 24:
        cond = tau_1rd1r(mat.shape[0], 0, bc.get('c',None))
    # D(1/r D(r u)) = 0
    elif bc[0] == 25:
        cond = tau_d1rdr(mat.shape[0], 0, bc.get('c',None))
    # u = Du = 0
    elif bc[0] == 40:
        cond = tau_value_diff(mat.shape[0], 0, bc.get('c',None))
    # u = D^2u = 0
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
        cond.append([c*norm_c(i) for i in np.arange(0,nr)])
        c = next(it)

    if pos <= 0:
        cond.append([c*norm_c(i)*(-1.0)**i for i in np.arange(0,nr)])

    return np.array(cond)

def tau_integral(nr, pos, coeffs = None):
    """Create the boundary integral tau line(s)"""

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
        tmp = np.zeros((1,nr))
        tmp[0,::2] = [2*(n/(n**2-1) - 1/(n-1)) for n in np.arange(0,nr,2)]
        tmp[0,0] = tmp[0,0]/2
        cond.append(tmp[0,:])
        c = next(it)

    if pos <= 0:
        tmp = np.zeros((1,nr))
        tmp[0,::2] = [2*(n/(n**2-1) - 1/(n-1)) for n in np.arange(0,nr,2)]
        tmp[0,0] = tmp[0,0]/2
        cond.append(tmp[0,:])

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

    return np.array(cond)

def tau_1rdr(nr, pos, coeffs = None):
    """Create the 1/r D(r u) tau line(s)"""

    a = coeffs['a']
    b = coeffs['b']
    cond = tau_diff(nr, pos, 1/a) + tau_value(nr, pos, [1/(a+b), 1/(a-b)])

    return cond

def tau_1rd1r(nr, pos, coeffs = None):
    """Create the 1/r D(u/r) tau line(s)"""

    a = coeffs['a']
    b = coeffs['b']
    cond = tau_diff(nr, pos, 1/a) - tau_value(nr, pos, [1/(a+b), 1/(a-b)])

    return cond

def tau_d1rdr(nr, pos, coeffs = None):
    """Create the D(1/r D(r u)) tau line(s)"""

    a = coeffs['a']
    b = coeffs['b']
    cond = tau_diff2(nr, pos, 1/a) + tau_diff(nr, pos, [1/(a+b), 1/(a-b)])

    return cond

def tau_value_diff(nr, pos, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    cond = []
    if pos >= 0:
        cond.append(list(tau_value(nr,1,coeffs)[0]))
        cond.append(list(tau_diff(nr,1,coeffs)[0]))

    if pos <= 0:
        cond.append(list(tau_value(nr,-1,coeffs)[0]))
        cond.append(list(tau_diff(nr,-1,coeffs)[0]))

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

    return np.array(cond)

def apply_galerkin(mat, bc):
    """Apply a Galerkin stencil on the matrix"""

    return mat

def norm_c(n):
    """Compute the chebyshev normalisation c factor"""

    if n > 0:
        return 2
    else:
        return 1
