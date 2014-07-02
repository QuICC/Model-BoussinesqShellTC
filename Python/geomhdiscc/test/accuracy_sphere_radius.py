"""Check accuracy for radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import geomhdiscc.geometry.cartesian.chebyshev_tools as ct
import geomhdiscc.geometry.spherical.sphere_radius as sphere

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
l_even = 2
l_odd = 3
ls = [l_even, l_odd]

#
# Accuracy tests: zblk
#
print("zblk:")
A = sphere.zblk(nr, 0, no_bc)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
ssol = 0
#test_forward(A, sphys, ssol, xg, 0)


#
# Accuracy tests: i2x2 even
#
print("i2x2:")
for l in ls:
    A = sphere.i2x2(nr, l, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
    sphys = sy.expand(sphys*x**2)
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 1)


#
# Accuracy tests: i2x2lapl
#
print("i2x2lapl:")
for l in ls:
    A = sphere.i2x2lapl(nr, l, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
    sphys = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1))
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 1)


#
# Accuracy tests: i4x4
#
print("i4x4:")
for l in ls:
    A = sphere.i4x4(nr, l, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
    sphys = sy.expand(sphys*x**4)
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)


#
# Accuracy tests: i4x4lapl
#
print("i4x4lapl:")
for l in ls:
    A = sphere.i4x4lapl(nr, l, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
    sphys = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2)
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)


#
# Accuracy tests: i4x4lapl2
#
print("i4x4lapl2:")
for l in ls:
    A = sphere.i4x4lapl2(nr, l, no_bc)
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
    sphys = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) - (l-1)*l*(l+1)*(l+2))
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 2)
