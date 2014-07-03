"""Check accuracy for radial direction in a sphere"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import geomhdiscc.transform.sphere as transf
import geomhdiscc.geometry.spherical.sphere_radius as sphere


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)


def test_forward(op, parity, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    try:
        pres, psol = parity
    except:
        pres = parity
        psol = parity

    x = sy.Symbol('x')
    lhs = transf.tocheb(x_to_phys(res_expr,grid), pres)
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.tocheb(t, psol)
    err = np.abs(rhs - sol)
    print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))


def zblk(nr, ls, bc, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, 0, no_bc)
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = 0.0
        test_forward(A, l%2, sphys, ssol, rg, 0)


def i2x2(nr, ls, bc, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i2x2(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)


def i2x2lapl(nr, ls, bc, rg):
    """Accuracy test for zblk operator"""

    print("i2x2lapl:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i2x2lapl(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)


def i4x4(nr, ls, bc, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i4x4(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)


def i4x4lapl(nr, ls, bc, rg):
    """Accuracy test for i4x4lapl operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i4x4lapl(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)


def i4x4lapl2(nr, ls, bc, rg):
    """Accuracy test for i4x4lapl2 operator"""

    print("i4x4lapl2:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i4x4lapl2(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)


def i2x1(nr, ls, bc, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i2x1(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)


def i2x2d1(nr, ls, bc, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i2x2d1(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)


def i4x3(nr, ls, bc, rg):
    """Accuracy test for i4x3 operator"""

    print("i4x3:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i4x3(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**3*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)


def i4x4d1(nr, ls, bc, rg):
    """Accuracy test for i4x4d1 operator"""

    print("i4x4d1:")
    x = sy.Symbol('x')
    for l in ls:
        print("\tTest for l = " + str(l))
        A = sphere.i4x4d1(nr, l, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)


def qid(nr, ls, bc, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    for l in ls:
        A = sphere.qid(nr, 3, bc)
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sphys
        test_forward(A, l%2, sphys, ssol, xg, 3)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    rg = transf.grid(nr)
    ls = [2, 3]
    no_bc = [0]

    # run tests
    #zblk(nr, ls, no_bc, rg)
    i2x2(nr, ls, no_bc, rg)
    i2x2lapl(nr, ls, no_bc, rg)
    i4x4(nr, ls, no_bc, rg)
    i4x4lapl(nr, ls, no_bc, rg)
    i4x4lapl2(nr, ls, no_bc, rg)
    i2x1(nr, ls, no_bc, rg)
    i2x2d1(nr, ls, no_bc, rg)
    i4x3(nr, ls, no_bc, rg)
    i4x4d1(nr, ls, no_bc, rg)
    qid(nr, ls, no_bc, rg)
