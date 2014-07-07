"""Module provides functions to generate the boundary conditions in a cartesian 1D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils


def constrain(mat, bc, eq_zrows):
    """Contrain the matrix with the (Tau or Galerkin) boundary condition"""

    if bc[0] > 0:
        bc_mat = apply_tau(mat, bc)
    elif bc[0] < 0:
        bc_mat = apply_galerkin(mat, bc, eq_zrows)
    else:
        bc_mat = mat

    return bc_mat


def apply_tau(mat, bc):
    """Add Tau lines to the matrix"""

    if bc[0] == 10:
        cond = tau_value(mat.shape[0], 1, bc[1:])
    elif bc[0] == 11:
        cond = tau_value(mat.shape[0], -1, bc[1:])
    elif bc[0] == 12:
        cond = tau_diff(mat.shape[0], 1, bc[1:])
    elif bc[0] == 13:
        cond = tau_diff(mat.shape[0], -1, bc[1:])
    elif bc[0] == 20:
        cond = tau_value(mat.shape[0], 0, bc[1:])
    elif bc[0] == 21:
        cond = tau_diff(mat.shape[0], 0, bc[1:])
    elif bc[0] == 40:
        cond = tau_value_diff(mat.shape[0], 0, bc[1:])
    elif bc[0] == 41:
        cond = tau_value_diff2(mat.shape[0], 0, bc[1:])

    if cond.dtype == 'complex_':
        bc_mat = mat.astype('complex_').tolil()
    else:
        bc_mat = mat.tolil()

    bc_mat[0:cond.shape[0],:] = cond

    return bc_mat


def tau_value(nx, pos, coeffs = None):
    """Create the tau line(s) for a zero boundary value"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    if pos >= 0:
        cond.append([c*tau_c(i) for i in np.arange(0,nx)])

    if pos <= 0:
        cond.append([c*tau_c(i)*(-1.0)**i for i in np.arange(0,nx)])

    return np.array(cond)


def tau_diff(nx, pos, coeffs = None):
    """Create the tau line(s) for a zero 1st derivative"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    if pos >= 0:
        cond.append([c*i**2 for i in np.arange(0,nx)])

    if pos <= 0:
        cond.append([(-1.0)**(i+1)*c*i**2 for i in np.arange(0,nx)])

    return np.array(cond)


def tau_diff2(nx, pos, coeffs = None):
    """Create the tau line(s) for a zero 2nd derivative"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    if pos >= 0:
        cond.append([c*((1/3)*(i**4 - i**2)) for i in np.arange(0,nx)])

    if pos <= 0:
        cond.append([c*(((-1.0)**i/3)*(i**4 - i**2)) for i in np.arange(0,nx)])

    return np.array(cond)


def tau_value_diff(nx, pos, coeffs = None):
    """Create the tau lines for a zero boundary value and a zero 1st derivative"""

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    if pos >= 0:
        cond.append(list(tau_value(nx,1,coeffs)[0]))
        cond.append(list(tau_diff(nx,1,coeffs)[0]))

    if pos <= 0:
        cond.append(list(tau_value(nx,-1,coeffs)[0]))
        cond.append(list(tau_diff(nx,-1,coeffs)[0]))

    return np.array(cond)


def tau_value_diff2(nx, pos, coeffs = None):
    """Create tau lines for a zero boundary value and a zero 2nd derivative """

    if coeffs is None or len(coeffs) < 1:
        c = 1.0
    else:
        c = coeffs[0]

    cond = []
    if pos >= 0:
        cond.append(list(tau_value(nx,1,coeffs)[0]))
        cond.append(list(tau_diff2(nx,1,coeffs)[0]))

    if pos <= 0:
        cond.append(list(tau_value(nx,-1,coeffs)[0]))
        cond.append(list(tau_diff2(nx,-1,coeffs)[0]))

    return np.array(cond)


def stencil(nx, bc):
    """Create a Galerkin stencil matrix"""

    if bc[0] == -10:
        mat = stencil_value(nx, 1)
    elif bc[0] == -11:
        mat = stencil_value(nx, -1)
    elif bc[0] == -12:
        mat = stencil_value(nx, 1)
    elif bc[0] == -13:
        mat = stencil_value(nx, -1)
    elif bc[0] == -20:
        mat = stencil_value(nx, 0)
    elif bc[0] == -21:
        mat = stencil_diff(nx, 0)
    elif bc[0] == -40:
        mat = stencil_value_diff(nx, 0)
    elif bc[0] == -41:
        mat = stencil_value_diff2(nx, 0)

    return mat


def apply_galerkin(mat, bc, eq_zero_rows):
    """Apply a Galerkin stencil on the matrix"""

    nx = mat.shape[0]

    return stencil_eye(nx, eq_zero_rows)*mat*stencil(nx, bc)


def stencil_eye(nx, q):
    """Create the restriction identity to resize matrix after stencil use"""

    offsets = [q]
    diags = [[1]*(nx-q)]

    return spsp.diags(diags, offsets, (nx-q, nx))


def stencil_value(nx, pos):
    """Create stencil matrix for a zero boundary value"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos 

    # Generate subdiagonal
    def d_1(n):
        return galerkin_c(n+offsets[0])*sgn

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets, (nx,nx+offsets[0]))


def stencil_diff(nx, pos):
    """Create stencil matrix for a zero 1st derivative"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]
        sgn = -1.0
    else:
        offsets = [-1, 0]
        sgn = -pos 

    # Generate subdiagonal
    def d_1(n):
        return sgn*(n+offsets[0])**2/n**2

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets)


def stencil_diff2(nx, pos):
    """Create stencil matrix for a zero 2nd derivative"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-2, 0]

        # Generate subdiagonal
        def d_1(n):
            return -(n - 3.0)*(n - 2.0)**2/(n**2*(n + 1.0))

    else:
        offsets = [-1, 0]

        # Generate subdiagonal
        def d_1(n):
            return -pos*(n - 2.0)*(n - 1.0)/(n*(n + 1.0))

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets)


def stencil_value_diff(nx, pos):
    """Create stencil matrix for a zero boundary value and a zero 1st derivative"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-4, -2, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            return (n - 3.0)/(n - 1.0)

        # Generate 1st subdiagonal
        def d_1(n):
            return -2.0*n/(n + 1.0)

    else:
        offsets = [-2, -1, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            return (2.0*n - 3.0)/(2.0*n - 1.0)

        # Generate 1st subdiagonal
        def d_1(n):
            return -pos*4.0*n/(2.0*n + 1.0)

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_2, d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets)


def stencil_value_diff2(nx, pos):
    """Create stencil matrix for a zero boundary value and a zero 2nd derivative"""

    ns = np.arange(0,nx,1)
    if pos == 0:
        offsets = [-4, -2, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            return (n - 3.0)*(2.0*n**2 - 12.0*n + 19.0)/((n - 1.0)*(2*n**2 - 4.0*n + 3.0))

        # Generate 1st subdiagonal
        def d_1(n):
            return -2.0*n*(2.0*n**2 + 7.0)/((n + 1.0)*(2.0*n**2 + 4.0*n + 3.0))

    else:
        offsets = [-2, -1, 0]

        # Generate 2nd subdiagonal
        def d_2(n):
            return (n - 3.0)*(2.0*n - 3.0)/(n*(2.0*n - 1.0))

        # Generate 1st subdiagonal
        def d_1(n):
            return -pos*2.0*(2*.0*n**2 + 1.0)/((n + 1.0)*(2.0*n + 1.0))

    # Generate diagonal
    def d0(n):
        return 1.0

    ds = [d_2, d_1, d0]
    diags = utils.build_diagonals(ns, -1, ds, offsets, None, False)
    diags[-1] = diags[-1][0:nx+offsets[0]]

    return spsp.diags(diags, offsets)


def tau_c(n):
    """Compute the chebyshev normalisation c factor for tau boundary"""

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
