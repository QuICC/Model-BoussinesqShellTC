# Nicolò Lardelli 24th November 2016

from __future__ import division
from __future__ import unicode_literals

from scipy import fftpack as fftpack
import numpy as np
from numpy.polynomial import legendre as leg
from numpy.polynomial import chebyshev as cheb
from scipy import special as spe
from projection.spherical import plm

from geomhdiscc.transform.spherical import thgrid, totphys, totleg, eqgrid

#min_r_points = 2000

"""
def rgrid(nr, a, b):
    #Create the radial Chebyshev grid

    gN = max(min_r_points, nr)

    return a * np.cos(np.pi * (np.arange(0, gN) + 0.5) / gN) + b

"""
"""
grid_1d = rgrid
grid_fast = rgrid
grid_slow = thgrid
"""

def grid_2d(nr, a, b, maxl, m):
    """Compute the 2D grid for the contours"""

    r = rgrid(nr, a, b)
    th = thgrid(maxl, m)
    rmesh, thmesh = np.meshgrid(r, th)
    X = rmesh * np.sin(thmesh)
    Y = rmesh * np.cos(thmesh)

    return (X, Y)


def grid_fast_per(nr, a, b, m):
    """Compute the 2D grid for the equatorial contours"""

    r = rgrid(nr, a, b)
    phi = eqgrid(m)
    rmesh, phimesh = np.meshgrid(r, phi)
    X = rmesh * np.cos(phimesh)
    Y = rmesh * np.sin(phimesh)

    return (X, Y)


def grid_slow_per(maxl, tm, m):
    """Compute the 2D grid for the equatorial contours"""

    th = thgrid(maxl, tm)
    phi = eqgrid(m)
    X, Y = np.meshgrid(th, phi)

    return (X, Y)


def torphys(spec):
    """Transform radial spectral coefficients to physical values"""

    if len(spec) < min_r_points:
        data = np.hstack((spec, np.zeros(min_r_points - len(spec))))
    else:
        data = spec

    return fftpack.dct(data, 3)


toprofile = torphys


def torcheb(phys):
    """Transform radial physical values to spectral coefficients"""

    n = len(phys)

    return fftpack.dct(phys, 2) / (2 * n)


def torphys2D(spec):
    """Transform 2D R spectral coefficients to 2D R physical values"""

    phys = np.zeros((max(spec.shape[0], min_r_points), spec.shape[1]))
    for j in range(spec.shape[1]):
        phys[:, j] = torphys(spec[:, j])

    return phys


def torcheb2D(phys):
    """Transform 2D R physical values to 2D R spectral coefficients"""

    for j in range(phys.shape[1]):
        phys[:, j] = torcheb(phys[:, j])

    return phys


def toslice(spec, nr, a, b, maxl, m):
    """Transform to latitudinal slice"""

    spec = np.reshape(spec, (nr, maxl - m + 1), order='F')

    rphys = torphys2D(spec.real) + 1j * torphys2D(spec.imag)

    rspec = rphys.T
    phys = totphys(rspec, maxl, m)
    return phys

def projphyslmn(l, m, n, x, eq_params):
    # generate the projection for the lmn basis to the points x (rows: r, theta, phi; columns: colums: different points)

    try:
        r = x[:,0]
        theta = x[:,1]
        phi = x[:,2]
        coeff = np.zeros(np.max(n)+1)
        coeff[n]=1.
    except RuntimeError as e:
        print(e)
        raise

    # retrieve the correct radial transform from eq_params or somewhere similar
    xi = eq_params['rratio']

    # transform the radius in the chebyshev domain; we consider r \in [xi,1]
    r_norm = 2./(1.-xi)*r -(1.+xi)/(1.-xi)
    x_cos = np.cos(theta)

    T = cheb.chebval(r_norm, coeff)

    P = plm(l, m, x_cos)

    Q = np.exp(1j*phi*m)
    # Q_s = np.sin(phi*m)

    return T*P*Q

def proj_radial(nr, a, b, x):
    # evaluate the radial basis functions of degree <nr

    coeffs = np.eye(nr)

    # TODO: perhaps do an assertion to check that the projected value is within [-1,1]
    return cheb.chebval(a*x+b,coeffs)/x**2

def proj_dradial_dr(nr,a,b,x):
     #evaluate the first derivatice of radial basis function of degree <nr

     coeffs = np.eye(nr)

     c = cheb.chebder(coeffs)
     return  cheb.chebval(a*x+b,c)*a/x




if __name__=="__main__":

    eq_params = {'rratio':0.333333}
    x = np.array([[0.8, 0.2, 0.0], [0.9, 0.2, 0.3]])
    #print(projphyslmn(10,2,4,x,eq_params))
    a = 3.
    b = -2.
    #print(proj_radial(20,a,b,x[:,0]))

    print(proj_dradial_dr(20,a,b,x[:,0]))

    print(proj_radial(20,a,b,x[:,0]).shape)

    print(proj_dradial_dr(20,a,b,x[:,0]).shape)

