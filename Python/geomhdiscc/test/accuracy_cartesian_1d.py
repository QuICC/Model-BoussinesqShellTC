"""Check accuracy for cartesian 1D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import geomhdiscc.transform.cartesian as transf
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_generic_1d as cg1d


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)


def test_forward(op, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.tocheb(x_to_phys(res_expr,grid))
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.tocheb(t)
    err = np.abs(rhs - sol)
    print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))


def zblk(nx, bc, xg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    A = c1d.zblk(nx)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = 0
    test_forward(A, sphys, ssol, xg, 0)


def d1(nx, bc, xg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    A = c1d.d1(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sy.diff(sphys,x)
    test_forward(A, sphys, ssol, xg, 1)


def i1(nx, bc, xg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    A = c1d.i1(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, xg, 1)


def i2(nx, bc, xg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    A = c1d.i2(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 2)


def i2d2(nx, bc, xg):
    """Accuracy test for i2d2 operator"""

    print("i2d2:")
    x = sy.Symbol('x')
    A = c1d.i2d2(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)


def i2d1(nx, bc, xg):
    """Accuracy test for i2d1 operator"""

    print("i2d1:")
    x = sy.Symbol('x')
    A = c1d.i2d1(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.integrate(sy.diff(sphys,x),x,x)
    test_forward(A, sphys, ssol, xg, 2)


def i2lapl(nx, bc, xg):
    """Accuracy test for i2lapl operator"""

    print("i2lapl:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i2lapl(nx, k, l, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys - l**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)


def i2laplh(nx, bc, xg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i2laplh(nx, k, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)


def i4(nx, bc, xg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    A = c1d.i4(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4d2(nx, bc, xg):
    """Accuracy test for i4d2 operator"""

    print("i4d2:")
    x = sy.Symbol('x')
    A = c1d.i4d2(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4d4(nx, bc, xg):
    """Accuracy test for i4d4 operator"""

    print("i4d4:")
    x = sy.Symbol('x')
    A = c1d.i4d4(nx, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4lapl(nx, bc, xg):
    """Accuracy test for i4lapl operator"""

    print("i4lapl:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl(nx, k, l, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys - l**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4laplh(nx, bc, xg):
    """Accuracy test for i4laplh operator"""

    print("i4laplh:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i4laplh(nx, k, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4lapl2(nx, bc, xg):
    """Accuracy test for i4lapl2 operator"""

    print("i4lapl2:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl2(nx, k, l, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) + k**4*sphys + l**4*sphys - 2*k**2*sy.diff(sphys,x,x) - 2*l**2*sy.diff(sphys,x,x) + 2*k**2*l**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def i4lapl2h(nx, bc, xg):
    """Accuracy test for i4lapl2h operator"""

    print("i4lapl2h:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i4lapl2h(nx, k, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) + k**4*sphys - 2*k**2*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)


def qid(nx, bc, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    A = c1d.qid(nx, 3, bc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)


def genxp(nx, q, p, bc, xg, ntrunc = -1):
    """Accuracy test for genxp operator"""

    print("genx"+ str(p) + ":")
    x = sy.Symbol('x')
    expr = 'x**' + str(p)
    A = cg1d.generic(nx, q, expr, x, bc, 1.0, ntrunc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sphys*x**p
    test_forward(A, sphys, ssol, xg, 1)


def genlinxp(nx, q, p, bc, xg, ntrunc = -1):
    """Accuracy test for gen_linxp operator"""

    print("genlinx"+ str(p) + ":")
    x = sy.Symbol('x')
    expr = '(1+x)**' + str(p)
    print(expr)
    A = cg1d.generic(nx, q, expr, x, bc, 1.0, ntrunc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sphys*(1+x)**p
    test_forward(A, sphys, ssol, xg, 1)


if __name__ == "__main__":
    # Set test parameters
    nx = 50
    xg = transf.grid(nx)
    no_bc = [0]

    # run tests
    print('Hard coded exact operators')
    #zblk(nx, no_bc, xg)
    d1(nx, no_bc, xg)
    i1(nx, no_bc, xg)
    i2(nx, no_bc, xg)
    i2d1(nx, no_bc, xg)
    i2d2(nx, no_bc, xg)
    i2lapl(nx, no_bc, xg)
    i2laplh(nx, no_bc, xg)
    i4(nx, no_bc, xg)
    i4d2(nx, no_bc, xg)
    i4d4(nx, no_bc, xg)
    i4lapl(nx, no_bc, xg)
    i4laplh(nx, no_bc, xg)
    i4lapl2(nx, no_bc, xg)
    i4lapl2h(nx, no_bc, xg)
    qid(nx, no_bc, xg)

    # run generic operator tests
    print('Generic operator generator')
    genxp(nx, 0, 1, no_bc, xg)
    genxp(nx, 0, 2, no_bc, xg)
    genxp(nx, 0, 3, no_bc, xg)
    genxp(nx, 0, 4, no_bc, xg)
    genxp(nx, 0, 9, no_bc, xg)
    genlinxp(nx, 0, 1, no_bc, xg)
    genlinxp(nx, 0, 2, no_bc, xg)
    genlinxp(nx, 0, 1.3, no_bc, xg)
    genlinxp(nx, 0, 1.3, no_bc, xg, 10)
