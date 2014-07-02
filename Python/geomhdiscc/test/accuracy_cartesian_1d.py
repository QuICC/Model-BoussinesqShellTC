"""Check accuracy of operations for cartesian 1D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy
import geomhdiscc.geometry.cartesian.chebyshev_tools as ct
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import scipy.sparse as spsp

x = sympy.Symbol('x')

def x_to_phys(expr, grid):
    func = sympy.utilities.lambdify(x, expr)
    return func(grid)


def test_forward(op, res_expr, sol_expr, grid, q):
    lhs = ct.tocheb(x_to_phys(res_expr,grid))
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = ct.tocheb(t)
    err = np.abs(rhs - sol)
    print(err)
    print(np.max(err[q:]))


nx = 20
xg = ct.grid(nx)

#
# Accuracy tests: zblk
#
print("zblk:")
A = c1d.zblk(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = 0
#test_forward(A, sphys, ssol, xg, 0)

#
# Accuracy tests: i2d2
#
print("i2d2:")
A = c1d.i2d2(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x),x,x)
test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i1
#
print("i1:")
A = c1d.i1(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
ssol = sympy.integrate(sphys,x)
test_forward(A, sphys, ssol, xg, 1)

#
# Accuracy tests: i2
#
print("i2:")
A = c1d.i2(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
ssol = sympy.integrate(sphys,x,x)
test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i2lapl
#
print("i2lapl:")
k, l = np.random.rand(2)*nx
A = c1d.i2lapl(nx, k, l)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x) - k**2*sphys - l**2*sphys,x,x)
test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i2laplh
#
print("i2laplh:")
k = np.random.ranf()*nx
A = c1d.i2laplh(nx, k)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x) - k**2*sphys,x,x)
test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i4
#
print("i4:")
A = c1d.i4(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)])
ssol = sympy.integrate(sphys,x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4d2
#
print("i4d2:")
A = c1d.i4d2(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x),x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4d4
#
print("i4d4:")
A = c1d.i4d4(nx)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x,x,x),x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4lapl
#
print("i4lapl:")
k, l = np.random.rand(2)*nx
A = c1d.i4lapl(nx, k, l)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x) - k**2*sphys - l**2*sphys,x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4laplh
#
print("i4laplh:")
k = np.random.ranf()*nx
A = c1d.i4laplh(nx, k)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x) - k**2*sphys,x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4lapl2
#
print("i4lapl2:")
k, l = np.random.rand(2)*nx
A = c1d.i4lapl2(nx, k, l)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x,x,x) + k**4*sphys + l**4*sphys - 2*k**2*sympy.diff(sphys,x,x) - 2*l**2*sympy.diff(sphys,x,x) + 2*k**2*l**2*sphys,x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: i4lapl2h
#
print("i4lapl2h:")
k = np.random.ranf()*nx
A = c1d.i4lapl2h(nx, k)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sympy.integrate(sympy.diff(sphys,x,x,x,x) + k**4*sphys - 2*k**2*sympy.diff(sphys,x,x),x,x,x,x)
test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: qid
#
print("qid:")
A = c1d.qid(nx, 3)
sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
ssol = sphys
test_forward(A, sphys, ssol, xg, 3)
