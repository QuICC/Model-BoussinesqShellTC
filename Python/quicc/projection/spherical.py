# -*- coding: utf-8 -*-
# Nicolò Lardelli 24th November 2016

from __future__ import division
from __future__ import unicode_literals

import numpy as np


def plm(l, m, x):

    """Compute the normalized associated legendre polynomial projection matrix"""
    if l<m or m<0 or l<0:
        return np.zeros_like(x)
    maxl = l
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

    # orthogonality as Schäffer 2013
    mat[:,0] = 1.0/np.sqrt(4.0*np.pi)
    sx = np.sqrt(1.0 - x**2)
    for i in range(1, m+1):
        mat[:,0] = -np.sqrt((2.0*i + 1.0)/(2.0*i))*sx*mat[:,0]

    return mat

def lplm(l, m, x):

    # return an orthonormal assoc legendre function
    # wapper for plm
    x = np.array(x)

    y = plm(l, m, x)
    return y

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

def dplm_1(l, m, x):
    # derivative associated legendre function
    # implementation Hollerbach style, stable on the poles
    # fully normalized

    x = np.array(x)
    if(l==0 and m==0):
        return np.zeros_like(x)
    y = -1./2*((l+m)*(l-m+1)*spe.lpmv(m-1,l,x)-spe.lpmv(m+1,l,x))

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5

def dplm_1(l, m, x):

    x = np.array(x)
    y = (l+1)*(l+m)*(spe.lpmv(m+1,l-2,x)+(l+m-2)*(l+m-1)*spe.lpmv(m-1,l-2,x))
    y -= l*(l-m+1)*(spe.lpmv(m+1,l,x)+(l+m)*(l+m+1)*spe.lpmv(m-1,l,x))

    return y * (spe.gamma(l - m + 1) / spe.gamma(l + m + 1) * (2 * l + 1) / 4 / np.pi) ** .5 /((2*l+1)*2*m)

def lplm_sin_1(l, m, x):


    x = np.array(x)
    if m!=0:
        y = -1/2/m*(spe.lpmv(m+1, l-1, x) + (l+m-1)*(l+m)*spe.lpmv(m-1, l-1,x))
    else:
        y = spe.lpmv(m, l, x)/(1-x**2)**.5

    return y*(spe.gamma(l-m+1)/spe.gamma(l+m+1)*(2*l+1)/4/np.pi)**.5


def deipm(l, m, phi):
    phi = np.array(phi)
    # compute the derivative  wrt phi of the azimuthal part of  Y_l^m
    return 1j*m*np.exp(m*phi*1j)

"""

def dplm(l, m, x):

    # returns the derivative associated legendre function
    # implementation Hollerbach style, stable on the poles
    # fully normalized
    x = np.array(x)
    y = -1./2*(((l+m)*(l-m+1))**0.5*plm(l,m-1,x)-((l-m)*(l+m+1))**.5*plm(l, m+1, x))

    return y


def lplm_sin(l, m, x):

    # return associated legendre function /sin_theta
    # implemented with the recurrence relation for P_l^m(x)/(1-x**2)**.5 when possible (aka m!=0)
    # fully normalized

    x = np.array(x)
    if m!=0:
        y = -1/2/m*(((l-m)*(l-m-1))**.5 * plm(l-1, m+1, x) + ((l+m)*(l+m-1))**.5 *plm(l-1, m-1, x))*((2*l+1)/(2*l-1))**.5
    else:
        y = plm(l, m, x)/(1-x**2)**.5

    return y

def eipm(l, m, phi):

    # returns the azimuthal part of a Ylm spherical harmonic
    # the spherical harmonic is fully normalized, but the normalization
    # is hidden in the associated legendre function part
    phi = np.array(phi)

    # compute the azimuthal part of spherical harmonics e^{imphi}
    return np.exp(m*phi*1j)


if __name__=="__main__":
    # main function, used for testing the routines

    theta=np.linspace(0,np.pi,100)
    x = np.cos(theta)
    phi = np.linspace(0, 2*np.pi,10)







