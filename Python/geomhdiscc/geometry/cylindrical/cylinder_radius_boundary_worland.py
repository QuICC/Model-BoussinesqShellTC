"""Module provides functions to generate the radial boundary conditions in a cylinder with Worland expansion in radius"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.worland.worland_basis as wb 


def no_bc():
    """Get a no boundary condition flag"""

    return {0:0}

def constrain(mat, m, bc, pad_zeros = 0, location = 't'):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, m, bc, pad_zeros = pad_zeros, location = location)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, m, bc)
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

def apply_tau(mat, m, bc, pad_zeros = 0, location = 't'):
    """Add Tau lines to the matrix"""
    
    nbc = bc[0]//10

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 11:
        cond = tau_diff(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 12:
        cond = tau_rdiffdivr(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 14:
        cond = tau_diff2(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 15:
        cond = tau_laplh(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 16:
        cond = tau_dlaplh(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 17:
        cond = tau_lapl2h(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 20:
        cond = tau_value_diff(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 21:
        cond = tau_value_diff2(mat.shape[0], m, bc.get('c',None))
    elif bc[0] == 22:
        cond = tau_value_laplh(mat.shape[0], m, bc.get('c',None))
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
    if pad_zeros > 0:
        cond = conc((np.zeros((pad_zeros,cond.shape[1])),cond))
    elif pad_zeros < 0:
        cond = conc((cond, np.zeros((pad_zeros,cond.shape[1]))))
    for i,c in enumerate(cond):
        mat.data = conc((mat.data, c))
        mat.row = conc((mat.row, [s+i]*mat.shape[1]))
        mat.col = conc((mat.col, np.arange(0,mat.shape[1])))

    return mat

def tau_value(nr, l, coeffs = None):
    """Create the boundary value tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_value(nr, l))

    return np.array(cond)

def tau_diff(nr, l, coeffs = None):
    """Create the 1st derivative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_diff(nr,l))

    return np.array(cond)

def tau_diff2(nr, m, coeffs = None):
    """Create the second deriviative tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_diff2(nr,m))

    return np.array(cond)

def tau_rdiffdivr(nr, m, coeffs = None):
    """Create the r D 1/r tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_rdiffdivr(nr,m))

    return np.array(cond)

def tau_lapl2h(nr, m, coeffs = None):
    """Create the no-slip tau line(s) for toroidal component on poloidal"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_lapl2h_cyl(nr,m))

    return np.array(cond)

def tau_value_diff(nr, m, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_value(nr,m))
    cond.append(c*wb.worland_diff(nr,m))

    return np.array(cond)

def tau_laplh(nr, m, coeffs = None):
    """Create the no-slip tau line(s) for toroidal component on toroidal"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_laplh_cyl(nr,m))

    return np.array(cond)

def tau_dlaplh(nr, m, coeffs = None):
    """Create the no-slip tau line(s) for toroidal component on toroidal"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_dlaplh_cyl(nr,m))

    return np.array(cond)

def tau_value_diff2(nr, m, coeffs = None):
    """Create the no penetration and no-slip tau line(s)"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_value(nr,m))
    cond.append(c*wb.worland_diff2(nr,m))

    return np.array(cond)

def tau_value_laplh(nr, m, coeffs = None):
    """Create the no penetration and no-slip tau line(s) for poloidal component"""

    if coeffs is None:
        c = 1.0
    else:
        c = coeffs

    cond = []
    cond.append(c*wb.worland_value(nr,m))
    cond.append(c*wb.worland_laplh_cyl(nr,m))

    return np.array(cond)

def tau_last(nr, nrow):
    """Create the last modes to zero tau line(s)"""

    cond = np.zeros((nrow, nr))
    for j in range(0, nrow):
        cond[j,nr-nrow+j] =1.0 

    return cond

def stencil(nr, m, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nr, m, bc.get('c',None))
    elif bc[0] == -11:
        mat = stencil_diff(nr, m, bc.get('c',None))
    elif bc[0] == -12:
        mat = stencil_rdiffdivr(nr, m, bc.get('c',None))
    elif bc[0] == -20:
        mat = stencil_value_diff(nr, m, bc.get('c',None))
    elif bc[0] == -21:
        mat = stencil_value_diff2(nr, m, bc.get('c',None))
    elif bc[0] < -1 and bc[0] > -5:
        mat = restrict_eye(nr, 'cr', -bc[0])

    return mat

def apply_galerkin(mat, m, bc):
    """Apply a Galerkin stencil on the matrix"""

    nr = mat.shape[0]
    mat = mat*stencil(nr, l, bc)
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

def stencil_value(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        return -2.0*n/(2.0*n - 1.0)

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_diff(nr, l, coeffs = None):
    """Create stencil matrix for a zero 1st derivative"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        num = -2.0*m*(4.0*(-1.0 + m)**2 + l*(-3.0 + 4.0*m))
        den = (-1.0 + 2.0*m)*(l + 4.0*l*m + 4.0*m**2)
        return num/den

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_rdiffdivr(nr, l, coeffs = None):
    """Create stencil matrix for a zero r D 1/r"""

    assert(coeffs is None)

    ns = np.arange(0,nr)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        num = -2.0*n*(3.0 - 3.0*l - 8.0*n + 4.0*l*n + 4.0*n**2)
        den = (-1.0 + 2.0*n)*(-1.0 + l + 4.0*l*n + 4.0*n**2)
        return num/den

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_insulating(nr, l, coeffs = None):
    """Create stencil matrix for a insulating boundary"""

    raise NotImplementedError("Boundary condition not yet implemented!")
    assert(coeffs.get('l', None) is not None)

    l = coeffs['l']

    ns = np.arange(parity,2*nr,2)
    offsets = [-1, 0]

    # Generate subdiagonal
    def d_1(n):
        val = -(l + n**2 - 4.0*n + 5.0)/(l + n**2 + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -(l + 1.0)/(l + 5.0)
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 1st derivative"""

    raise NotImplementedError("Boundary condition not yet implemented!")
    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-2, -1, 0]

    # Generate second subdiagonal
    def d_2(n):
        val = (n - 3.0)/(n - 1.0)
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/6.0
                break
            if j > 4:
                break
        return val

    # Generate subdiagonal
    def d_1(n):
        val = -2.0*n/(n + 1.0)
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -2.0/3.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))

def stencil_value_diff2(nr, l, coeffs = None):
    """Create stencil matrix for a zero boundary value and zero 2nd derivative"""

    raise NotImplementedError("Boundary condition not yet implemented!")
    assert(coeffs is None)

    ns = np.arange(parity,2*nr,2)
    offsets = [-2, -1, 0]

    # Generate second subdiagonal
    def d_2(n):
        val_num = (n - 3.0)*(2.0*n**2 - 12.0*n + 19.0)
        val_den = (n - 1.0)*(2.0*n**2 - 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 4:
                val[i] = 1.0/38.0
                break
            if j > 4:
                break
        return val

    # Generate subdiagonal
    def d_1(n):
        val_num = -2.0*n*(2.0*n**2 + 7.0)
        val_den = (n + 1.0)*(2.0*n**2 + 4.0*n + 3.0)
        val = val_num/val_den
        for i,j in enumerate(n):
            if j == 2:
                val[i] = -10.0/19.0
                break
            if j > 2:
                break
        return val

    # Generate diagonal
    def d0(n):
        return np.ones(n.shape)

    ds = [d_2, d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nr+offsets[0]]

    return spsp.diags(diags, offsets, (nr,nr+offsets[0]))
