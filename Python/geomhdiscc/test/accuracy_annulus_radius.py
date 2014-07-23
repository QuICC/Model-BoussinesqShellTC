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
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))

def zblk(nr, ms, a, b, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.zblk(nr, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
        ssol = 0
        test_forward(A, sphys, ssol, rg, 0)

def d1(nr, ms, a, b, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.d1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.diff(sphys,x)
        test_forward(A, sphys, ssol, rg, 0)

def i1(nr, ms, a, b, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.integrate(sphys,x)
        test_forward(A, sphys, ssol, rg, 1)

def i1x1d1(nr, ms, a, b, rg):
    """Accuracy test for i1x1d1 operator"""

    print("i1x1d1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i1x1d1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sy.diff(sphys,x)*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, sphys, ssol, rg, 1)

def i1x1(nr, ms, a, b, rg):
    """Accuracy test for i1x1 operator"""

    print("i1x1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i1x1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, sphys, ssol, rg, 1)

def i2(nr, ms, a, b, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.integrate(sphys,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x1(nr, ms, a, b, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x2d1(nr, ms, a, b, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2d1(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sy.diff(sphys,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x2d2(nr, ms, a, b, rg):
    """Accuracy test for i2x2d2 operator"""

    print("i2x2d2:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2d2(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sy.diff(sphys,x,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x2(nr, ms, a, b, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x2div(nr, ms, a, b, rg):
    """Accuracy test for i2x2div operator"""

    print("i2x2div:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2div(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x*sphys)
        ssol = sy.expand(x*sy.diff(ssol,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i2x2laplh(nr, ms, a, b, rg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i2x2laplh(nr, m, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, sphys, ssol, rg, 2)

def i4x4(nr, ms, a, b, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4(nr, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)

def i4x4laplh(nr, ms, a, b, rg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4laplh:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4laplh(nr, m, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl2h(nr, ms, a, b, rg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    x = sy.Symbol('x')
    for m in ms:
        print("\tTest for m = " + str(m))
        A = annulus.i4x4lapl2h(nr, m, a, b, annulus.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, sphys, ssol, rg, 4)

def qid(nr, a, b, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    A = annulus.qid(nr, 3, annulus.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    ms = [2, 3]
    a = 0.325
    b = 0.675
    rg = (a*transf.grid(nr) + b)

    # run tests
    #zblk(nr, ms, a, b, rg)
    d1(nr, ms, a, b, rg)
    i1(nr, ms, a, b, rg)
    i1x1d1(nr, ms, a, b, rg)
    i1x1(nr, ms, a, b, rg)
    i2(nr, ms, a, b, rg)
    i2x1(nr, ms, a, b, rg)
    i2x2(nr, ms, a, b, rg)
    i2x2d2(nr, ms, a, b, rg)
    i2x2d1(nr, ms, a, b, rg)
    i2x2div(nr, ms, a, b, rg)
    i2x2laplh(nr, ms, a, b, rg)
    i4x4(nr, ms, a, b, rg)
    i4x4laplh(nr, ms, a, b, rg)
    i4x4lapl2h(nr, ms, a, b, rg)
    qid(nr, a, b, rg)
