# Nicolò Lardelli 24th November 2016

from __future__ import division
from __future__ import unicode_literals
import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial import legendre as leg
from quicc.projection.spherical import plm
from quicc.transform.spherical import thgrid, totphys, totleg, eqgrid
from scipy import fftpack as fftpack
from scipy import special as spe


#min_r_points = 2000



def proj_radial(nr, a, b, x):
    # evaluate the radial basis functions of degree <nr
    xx = np.array(x)

    # use chebyshev polynomial normalized for fft
    coeffs = np.eye(nr)*2
    coeffs[0,0]=1.
    print(cheb.chebval(xx,coeffs))

    #print(xx*a+b)
    # TODO: perhaps do an assertion to check that the projected value is within [-1,1]
    return (cheb.chebval(xx,coeffs)/(a*xx+b)**2).transpose()

def proj_dradial_dr(nr,a,b,x):
     #evaluate the first derivatice of radial basis function of degree <nr
     xx = np.array(x)
     coeffs = np.eye(nr) * 2
     coeffs[0, 0] = 1.
     c = cheb.chebder(coeffs)
     return  (cheb.chebval(xx,c)/a /(a*xx+b)).transpose()



def proj_radial_r(nr, a, b, x):
    # evaluate the radial basis functions of degree <nr
    xx = np.array(x)
    coeffs = np.eye(nr)*2
    coeffs[0,0]=1.

    temp = (cheb.chebval(xx,coeffs)/(a*xx+b)).transpose()
    return temp



if __name__=="__main__":

    eq_params = {'rratio':0.333333}
    #x = np.array([[0.8, 0.2, 0.0], [0.9, 0.2, 0.3]])
    x1 = np.array([0.9,0.95,1.])
    #x2 = np.array([0.7,0.8,0.9,1.])
    x2 =np.array([-1,1])
    #print(projphyslmn(10,2,4,x,eq_params))
    (a,b) = (0.5,1.0)
    #print(proj_radial(20,a,b,x[:,0]))

    #print(proj_dradial_dr(20,a,b,x[:,0]))

    #print(proj_radial(80,a,b,x[:,0]))

    #print(proj_radial(10,a,b,x1))
    tupla = (10,a,b,x2)
    print(proj_radial(10,a,b,x2))
    print(proj_dradial_dr(tupla))
    print(proj_radial_r(tupla))

    #print(proj_dradial_dr(20,a,b,x[:,0]).shape)

    #print(proj_dradial_dr(80,a,b,x[:,0]).dot(np.random.rand(80)))


