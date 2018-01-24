# -*- coding: utf-8 -*-
# Nicolo Lardelli 24th November 2016

from __future__ import division
from __future__ import unicode_literals
import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy import sparse as spsp

import quicc.base.utils as utils

def proj_radial(nr, a, b, x):

    # evaluate the radial basis functions of degree <nr
    # this projection operator includes the 1/r**2 part and corresponds
    # to the projection of the poloidal field in r direction (Q-part of a qst decomposition
    xx = np.array(x)

    # use chebyshev polynomials normalized for FCT
    coeffs = np.eye(nr)*2
    coeffs[0,0]=1.

    return (cheb.chebval(xx,coeffs)/(a*xx+b)**2).transpose()

def proj_dradial_dr(nr,a,b,x):

     # evaluate the first derivative of radial basis function of degree <nr
     # this projection operator corresponds to 1/r d/dr projecting from spectral to physical
     # used for the theta and phi contribution of the poloidal field (S-part of a qst decomposition)
     xx = np.array(x)

     # use chebyshev polynomials normalized for FCT
     coeffs = np.eye(nr) * 2
     coeffs[0, 0] = 1.
     c = cheb.chebder(coeffs)
     return  (cheb.chebval(xx,c)/a /(a*xx+b)).transpose()



def proj_radial_r(nr, a, b, x):

    # evaluate the radial basis functions of degree <nr over r
    # this projection operator includes the 1/r part and corresponds
    # to the a projection from spectral to physical for the toroidal field (T-part of a qst decomposition)
    xx = np.array(x)

    # use chebyshev polynomials normalized for FCT
    coeffs = np.eye(nr)*2
    coeffs[0,0]=1.

    temp = (cheb.chebval(xx,coeffs)/(a*xx+b)).transpose()
    return temp

def proj(nr, a, b, r):

    # evaluates the projection matrix matrix for the chebyshev basis
    xx = (np.array(r)-b)/a
    coeffs = np.eye(nr)*2
    coeffs[0,0]=1.

    return spsp.lil_matrix(cheb.chebval(xx, coeffs).transpose())


if __name__=="__main__":

    eq_params = {'rratio':0.333333}
    #x = np.array([[0.8, 0.2, 0.0], [0.9, 0.2, 0.3]])
    x1 = np.array([0.9,0.95,1.])
    #x2 = np.array([0.7,0.8,0.9,1.])
    x2 =np.array([-1,1])
    (a,b) = (0.5,1.0)




