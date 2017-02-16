# Nicolò Lardelli 24th November 2016

from __future__ import division
from __future__ import unicode_literals

import numpy as np

from numpy.polynomial import legendre as leg
import scipy.special as spe


def totphys(spec, maxl, m, x):
    """Tranform theta spectral coefficients to physical values"""

    mat = plm(maxl, m, x)
    phys = mat.dot(spec)

    return phys

def totleg(phys, maxl, m):
    """Tranform theta physical values to spectral coefficients"""

    nt = 3*(maxl - m + 1)//2
    x, w = leg.leggauss(nt)
    mat = plm(maxl, m).T
    spec = mat.dot(np.diag(w).dot(phys))

    return spec

def eqgrid(m, phi = 2*np.pi):
    """Create a equatorial (phi) grid for given harmonic order"""

    return np.linspace(0, phi, max(min_phi_points,3*m))

def plm(maxl, m, x):
    """Compute the normalized associated legendre polynomial projection matrix"""

    mat = np.zeros((len(x), maxl - m + 1))
    mat[:,0] = pmm(maxl, m, x)[:,0]
    if maxl == m:
        return mat[:,-1]

    mat[:,1] = np.sqrt(2.0*m + 3.0)*x*mat[:,0]
    for i, l in enumerate(range(m+2, maxl + 1)):
        mat[:,i+2] = np.sqrt((2.0*l + 1)/(l - m))*np.sqrt((2.0*l - 1.0)/(l + m))*x*mat[:,i+1] - np.sqrt((2.0*l + 1)/(2.0*l - 3.0))*np.sqrt((l + m - 1.0)/(l + m))*np.sqrt((l - m - 1.0)/(l - m))*mat[:,i]

    return mat[:,-1]

def pmm(maxl, m, x):
    """Compute the normalized associated legendre polynomial of order and degree m"""

    mat = np.zeros((len(x), 1))
    mat[:,0] = 1.0/np.sqrt(2.0)
    sx = np.sqrt(1.0 - x**2)
    for i in range(1, m+1):
        mat[:,0] = -np.sqrt((2.0*i + 1.0)/(2.0*i))*sx*mat[:,0]

    return mat

def lplm(l, m, x):
    # return an orthonormal assoc legendre func

    x = np.array(x)
    y = spe.lpmv(m, l, x)

    return y*(spe.gamma(l-m+1)/spe.gamma(l+m+1)*(2*l+1)/4/np.pi)**.5

"""
def dplm0(l, m, x):
    # return the deriavative of an orthonormal assoc legendre func
    # implementation Hollerbach style
    sin_the = (1 - x ** 2) ** .5

    y = spe.lpmv(m,l+1,x)*l*(l+1-m) -spe.lpmv(m,l-1,x)*(l+1)*(l+m)/(2*l+1)

    return y*(spe.gamma(l-m+1)/spe.gamma(l+m+1)*(2*l+1)/4/np.pi)**.5/sin_the/(2*l+1)

def dplm1(l, m, x):

    sin_the = (1 - x ** 2) ** .5

    y = l* x * spe.lpmv(m,l,x) - (l+m)*spe.lpmv(m,l-1,x)

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /sin_the

def dplm2(l, m, x):

    sin_the = (1 - x ** 2) ** .5

    y = -(l+1) * x * spe.lpmv(m,l,x) +(l-m+1)*spe.lpmv(m,l+1,x)

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /sin_the

def dplm3(l, m, x):

    sin_the = (1 - x ** 2) ** .5

    y = sin_the* spe.lpmv(m+1,l,x) +m*x*spe.lpmv(m,l,x)

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /sin_the

def dplm4(l, m, x):

    sin_the = (1 - x ** 2) ** .5

    y = -(l+m) * (l-m+1)*sin_the * spe.lpmv(m-1,l,x) -m *x*spe.lpmv(m,l,x)

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /sin_the
"""

def dplm(l, m, x):
    # derivative associated legendre function
    # implementation Hollerbach style, stable on the poles
    # fully normalized

    x = np.array(x)
    if(l==0 and m==0):
        return np.zeros_like(x)
    y = -1./2*((l+m)*(l-m+1)*spe.lpmv(m-1,l,x)-spe.lpmv(m+1,l,x))

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5

"""
def dplm_1(l, m, x):


    x = np.array(x)
    y = (l+1)*(l+m)*(spe.lpmv(m+1,l-2,x)+(l+m-2)*(l+m-1)*spe.lpmv(m-1,l-2,x))

    y -= l*(l-m+1)*(spe.lpmv(m+1,l,x)+(l+m)*(l+m+1)*spe.lpmv(m-1,l,x))
    #print(l,m,'dplm')

    #print(y)

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /((2*l+1)*2*m)
"""
def lplm_sin(l, m, x):


    x = np.array(x)
    if m!=0:
        y = -1/2/m*(spe.lpmv(m+1, l-1, x) + (l+m-1)*(l+m)*spe.lpmv(m-1, l-1,x))
    else:
        y = spe.lpmv(m, l, x)/(1-x**2)**.5

    #print(l,m,'lplm_sin')
    return y*(spe.gamma(l-m+1)/spe.gamma(l+m+1)*(2*l+1)/4/np.pi)**.5


def eipm(l, m, phi):

    # convert phi argument into a numpy array
    phi = np.array(phi)

    # compute the azimuthal part of spherical harmonics e^{imphi}
    return np.exp(m*phi*1j)

def deipm(l, m, phi):
    phi = np.array(phi)
    # compute the derivative  wrt phi of the azimuthal part of  Y_l^m
    return 1j*m*np.exp(m*phi*1j)

if __name__=="__main__":
    # main function, used for testing the routines

    theta=np.linspace(0,np.pi,100)
    x = np.cos(theta)
    phi = np.linspace(0, 2*np.pi,10)
    #print(lplm(1, 1, x))

    print(dplm(1,0,x))
    #print(dplm(3,0,x))
    #print(dplm(7,0,x))
    #print(lplm_sin(3,0,x))
    #print(lplm(2,0,x)/(1-x**2)**.5)
    #print(lplm_sin(2,0,x))
    #print(lplm(2, 0, x) / (1 - x ** 2) ** .5 -lplm_sin(2,0,x))
    #print(dplm(1, 1, x))
    #print(eipm(1, 1, phi))
    #print(deipm(1, 1, phi))

    """
    print(dplm1(1, 1, x))
    print(dplm2(1, 1, x))
    print(dplm3(1, 1, x))
    print(dplm4(1, 1, x))
    print(dplm5(1, 1, x))

    print(dplm(1,1,x)-dplm5(1,1,x))
    """




