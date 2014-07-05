"""Module provides functions to generate generic sparse operators for a single cartesian direction"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.linalg as splin
import scipy.sparse as spsp
import geomhdiscc.base.utils as utils
import geomhdiscc.transform.cartesian as transf
import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc


def generic(nx, q, expr, var, bc, coeff = 1.0):
    """Compute the spectral operator for a generic expression"""

    # Integrate expression q times
    #if q > 0:
    #    nvar = (var,)
    #    for i in range(1,q):
    #        nvar = nvar + (var,)
    #    expr = sy.integrate(expr, nvar)
    func = sy.utilities.lambdify(var, expr)

    # Convert to physical space values
    grid = transf.grid(2*nx)
    phys = func(grid)
    spec = transf.tocheb(phys)[0:nx]
    spec *= 2*(np.abs(spec) > np.spacing(1))
    spec[0] = spec[0]/2

    # Form operator matrix
    h = splin.hankel(spec)
    h[0,:] = 0
    mat = splin.toeplitz(np.append(2*spec[0], spec[1:]))
    mat = 0.5*(mat + h)

    mat = coeff*spsp.coo_matrix(mat)
    return c1dbc.constrain(mat, bc, 2)
