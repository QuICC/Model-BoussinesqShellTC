"""Check accuracy for radial direction in a cylindrical annulus"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import geomhdiscc.transform.annulus as transf
import geomhdiscc.geometry.cylindrical.annulus_radius as annulus


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


def zblk(nr, ms, a, b, bc, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.zblk(nr, 0, no_bc)
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        ssol = 0
        test_forward(A, sphys, ssol, rg, 0)


def i2x2(nr, ms, a, b, bc, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2(nr, a, b, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)


def i2x2lapl(nr, ms, a, b, bc, rg):
    """Accuracy test for zblk operator"""

    print("i2x2lapl:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2lapl(nr, m, a, b, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)


def i4x4(nr, ms, a, b, bc, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4(nr, a, b, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)


def i4x4lapl(nr, ms, a, b, bc, rg):
    """Accuracy test for i4x4lapl operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4lapl(nr, m, a, b, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)


def i4x4lapl2(nr, ms, a, b, bc, rg):
    """Accuracy test for i4x4lapl2 operator"""

    print("i4x4lapl2:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4lapl2(nr, m, a, b, no_bc)
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    ms = [2, 3]
    no_bc = [0]
    a = 0.325
    b = 0.675
    rg = (a*transf.grid(nr) + b)

    # run tests
    #zblk(nr, ms, a, b, no_bc, rg)
    i2x2(nr, ms, a, b, no_bc, rg)
    i2x2lapl(nr, ms, a, b, no_bc, rg)
    i4x4(nr, ms, a, b, no_bc, rg)
    i4x4lapl(nr, ms, a, b, no_bc, rg)
    i4x4lapl2(nr, ms, a, b, no_bc, rg)
