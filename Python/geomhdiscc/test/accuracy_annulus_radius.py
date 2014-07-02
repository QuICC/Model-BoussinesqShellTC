"""Check accuracy for radial direction in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import geomhdiscc.geometry.cartesian.chebyshev_tools as ct
import geomhdiscc.geometry.cylindrical.annulus_radius as annulus

x = sy.Symbol('x')

def x_to_phys(expr, grid):
    func = sy.utilities.lambdify(x, expr)
    return func(grid)


def test_forward(op, res_expr, sol_expr, grid, q):
    lhs = ct.tocheb(x_to_phys(res_expr,grid))
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = ct.tocheb(t)
    err = np.abs(rhs - sol)
    print(err)
    print(np.max(err[q:]))


nr = 20
xg = ct.grid(nr)
no_bc = [0]
ms = [2, 3]
a = 0.325
b = 0.675

#
# Accuracy tests: zblk
#
print("zblk:")
A = annulus.zblk(nr, 0, no_bc)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
ssol = 0
#test_forward(A, sphys, ssol, xg, 0)


#
# Accuracy tests: i2x2 even
#
print("i2x2:")
for m in ms:
    A = annulus.i2x2(nr, a, b, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    sphys = sy.expand(sphys*x**2)
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 1)


#
# Accuracy tests: i2x2lapl
#
print("i2x2lapl:")
for m in ms:
    A = annulus.i2x2lapl(nr, m, a, b, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    sphys = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2)
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 1)


#
# Accuracy tests: i4x4
#
print("i4x4:")
for m in ms:
    A = annulus.i4x4(nr, a, b, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    sphys = sy.expand(sphys*x**4)
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)


#
# Accuracy tests: i4x4lapl
#
print("i4x4lapl:")
for m in ms:
    A = annulus.i4x4lapl(nr, m, a, b, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    sphys = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2)
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)


#
# Accuracy tests: i4x4lapl2
#
print("i4x4lapl2:")
for m in ms:
    A = annulus.i4x4lapl2(nr, m, a, b, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    sphys = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) - (m**2 - 4)*m**2)
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)
