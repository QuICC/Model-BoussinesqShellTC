"""Module provides functions to generate values for the Worland expansion"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import scipy.special as special

def worland_norm(n , l):
    """Normalization factor"""

    if l == 0 and n == 0:
        return np.sqrt(np.pi/2.0)
    else:
        return np.sqrt((2.0*n+l)*np.exp(special.gammaln(n+l) + special.gammaln(n+1.0) - special.gammaln(n+0.5) - special.gammaln(n+l+0.5)))

def worland_grid(nr):
    """Physical space grid"""

    return np.sqrt((np.cos(np.pi*(np.arange(0,2*nr)+0.5)/(2*nr)) + 1.0)/2.0)

def worland_weight(i, nr):
    """Gaussian integration weight"""

    return np.pi/(2.0*nr)

def worland_poly(n, l, nr):

    rg = worland_grid(nr)
    norm = worland_norm(n, l)

    return rg**l*special.eval_jacobi(n, -0.5, l - 0.5, 2.0*rg**2 - 1.0)*norm

def worland_norm_row(n , l, k):
    """Normalization factor for matrix row"""
    
    norm = np.log(2.0*n + l + 2.0*k) - np.log(2.0*n + l)
    norm += special.gammaln(n + 0.5) - special.gammaln(n + 0.5 + k)
    norm += special.gammaln(n + l + 0.5) - special.gammaln(n + l + 0.5 + k)
    norm += special.gammaln(n + l + k) - special.gammaln(n + l)
    norm += special.gammaln(n + 1.0 + k) - special.gammaln(n + 1.0)

    return np.sqrt(np.exp(norm))

def worland_normalize(val, l):
    """Normalize array of values"""

    for i in range(0, val.shape[0]):
         val[i] *= worland_norm(i, l)

def worland_value(nr, l, k = 0, normalize = True):
    """Compute the endpoint value for Worland polynomials"""

    val = np.zeros(nr)
    val[0] = 1.0
    if nr > 0:
        for i in range(1,nr):
            val[i] = val[i-1]*(2.0*i-1.0 + 2*k)/(2.0*i)

    # Normalize
    if normalize:
        worland_normalize(val, l)

    return val

def worland_diff(nr, l):
    """Compute the first derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[1:] = 2.0*(l+np.arange(1,nr))*worland_value(nr-1, l+1, 1, False)
        if l > 0:
            val += l*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_diff2(nr, l):
    """Compute the second derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[2:] = 4.0*(l+np.arange(2,nr)+1)*(l+np.arange(2,nr))*worland_value(nr-2, l+2, 2, False)
        if l > 0:
            val[1:] += (2.0*(l+np.arange(1,nr)) + 4.0*l*(l+np.arange(1,nr)))*worland_value(nr-1, l+1, 1, False)
        else:
            val[1:] += 2.0*(l+np.arange(1,nr))*worland_value(nr-1, l+1, 1, False)
        if l > 1:
            val += l*(l-1.0)*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_diff3(nr, l):
    """Compute the third derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[3:] = 8.0*(l+np.arange(3,nr)+2)*(l+np.arange(3,nr)+1)*(l+np.arange(3,nr))*worland_value(nr-3, l+3, 3, False)
        val[2:] += 12.0*(l+np.arange(2,nr))*(l+np.arange(2,nr)+1.0)*(1.0 + l)*worland_value(nr-2, l+2, 2, False)
        if l > 0:
            val[1:] += 6.0*l*(l+np.arange(1,nr))*l*worland_value(nr-1, l+1, 1, False)
        if l > 2:
            val += (l-1.0)*(l-2.0)*l*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_diff4(nr, l):
    """Compute the fourth derivative at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[4:] = 16.0*(l+np.arange(4,nr)+3)*(l+np.arange(4,nr)+2)*(l+np.arange(4,nr)+1)*(l+np.arange(4,nr))*worland_value(nr-4, l+4, 4, False)
        val[3:] += 16.0*(2.0*l + 3.0)*(l+np.arange(3,nr))*(l+np.arange(3,nr)+1.0)*(l+np.arange(3,nr)+2.0)*worland_value(nr-3, l+3, 3, False)
        val[2:] += 12.0*(2.0*l**2 + 2.0*l + 1.0)*(l+np.arange(2,nr))*(l+np.arange(2,nr)+1)*worland_value(nr-2, l+2, 2, False)
        if l > 1:
            val[1:] += 4.0*(l-1.0)*l*(2.0*l-1.0)*(l+np.arange(1,nr))*worland_value(nr-1, l, 1, False)
        if l > 3:
            val += l*(l-1.0)*(l-2.0)*(l-3.0)*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_rdiffdivr(nr, l):
    """Compute the stress-free condition at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[1:] = 2.0*(l+np.arange(1,nr))*worland_value(nr-1, l+1, 1, False)
        val += (l-1.0)*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_divrdiffr(nr, l):
    """Compute the insulating magnetic condition at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[1:] = 2.0*(l+np.arange(1,nr))*worland_value(nr-1, l+1, 1, False)
        val += (l+1.0)*worland_value(nr, l, 0, False)

    # Normalize
    worland_normalize(val, l)

    return val

def worland_laplh_cyl(nr, m):
    """Compute the horizontal laplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[2:] = 4.0*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*worland_value(nr-2, m+2, 2, False)
        val[1:] += 4.0*(m+1.0)*(m+np.arange(1,nr))*worland_value(nr-1, m+1, 1, False)

    # Normalize
    worland_normalize(val, m)

    return val

def worland_dlaplh_cyl(nr, m):
    """Compute the radial derivative of horizontal laplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[3:] = 8.0*(m+np.arange(3,nr)+2)*(m+np.arange(3,nr)+1)*(m+np.arange(3,nr))*worland_value(nr-3, m+3, 3, False)
        val[2:] += 4.0*(3.0*m + 4.0)*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*worland_value(nr-2, m+2, 2, False)
        if m > 0:
            val[1:] += 4.0*m*(m + 1.0)*(m+np.arange(1,nr))*worland_value(nr-1, m+1, 1, False)

    # Normalize
    worland_normalize(val, m)

    return val

def worland_lapl2h_cyl(nr, m):
    """Compute the bilaplacian in a cylinder at endpoint for Worland polynomials"""

    val = np.zeros(nr)
    if nr > 0:
        val[4:] = 16.0*(m+np.arange(4,nr)+3)*(m+np.arange(4,nr)+2)*(m+np.arange(4,nr)+1)*(m+np.arange(4,nr))*worland_value(nr-4, m+4, 4, False)
        val[3:] += 32.0*(m + 2.0)*(m+np.arange(3,nr)+2)*(m+np.arange(3,nr)+1)*(m+np.arange(3,nr))*worland_value(nr-3, m+3, 3, False)
        val[2:] += 16.0*(m + 1.0)*(m + 2.0)*(m+np.arange(2,nr)+1)*(m+np.arange(2,nr))*worland_value(nr-2, m+2, 2, False)

    # Normalize
    worland_normalize(val, m)

    return val
