"""Check accuracy for radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp

import geomhdiscc.transform.cylinder as transf
import geomhdiscc.geometry.cylindrical.cylinder_radius as cylinder


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
    lhs = transf.torcheb(x_to_phys(res_expr,grid), pres)
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torcheb(t, psol)
    err = np.abs(rhs - sol)
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))

def zblk(nr, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.zblk(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(parity,2*nr,2)])
        ssol = 0
        test_forward(A, parity, sphys, ssol, rg, 0)

def x1(nr, rg):
    """Accuracy test for x1 operator"""

    print("x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x*sphys)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def d1(nr, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x))
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def x1div(nr, rg):
    """Accuracy test for x1div operator"""

    print("x1div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.x1div(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(x*sphys,x))
        test_forward(A, parity, sphys, ssol, rg, 1)

def i1(nr, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.integrate(sphys,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1d1(nr, rg):
    """Accuracy test for i1x1d1 operator"""

    print("i1x1d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1div(nr, rg):
    """Accuracy test for i1x1div operator"""

    print("i1x1div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1div(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(x*sphys,x))
        ssol = sy.integrate(ssol,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i1x1(nr, rg):
    """Accuracy test for i1x1 operator"""

    print("i1x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i1x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i2(nr, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.integrate(sphys,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i2x1(nr, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i2x2d1(nr, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2d1(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 2)

def i2x2d2(nr, rg):
    """Accuracy test for i2x2d2 operator"""

    print("i2x2d2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2d2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sy.diff(sphys,x,x)*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i2x2(nr, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i2x2div(nr, rg):
    """Accuracy test for i2x2div operator"""

    print("i2x2div:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2div(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.expand(x*sy.diff(ssol,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (parity,(parity+1)%2), sphys, ssol, rg, 1)

def i2x2laplh(nr, rg):
    """Accuracy test for i2x2laplh operator"""

    print("i2x2laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i2x2laplh(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + x*sy.diff(sphys,x) - m**2*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, parity, sphys, ssol, rg, 1)

def i4(nr, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4(nr, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4(nr, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4laplh(nr, rg):
    """Accuracy test for i4x4laplh operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4laplh(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def i4x4lapl2h(nr, rg):
    """Accuracy test for i4x4lapl2h operator"""

    print("i4x4lapl2h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        print("\tTest for m = " + str(m))
        A = cylinder.i4x4lapl2h(nr, m, parity, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(parity,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, parity, sphys, ssol, rg, 2)

def qid(nr, rg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        parity = m%2
        A = cylinder.qid(nr, parity, 3, cylinder.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(m%2,2*nr,2)])
        ssol = sphys
        test_forward(A, parity, sphys, ssol, rg, 3)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    rg = transf.rgrid(nr)

    # run tests
    #zblk(nr, rg)
    x1(nr, rg)
    d1(nr, rg)
    x1div(nr, rg)
    i1(nr, rg)
    i1x1d1(nr, rg)
    i1x1div(nr, rg)
    i1x1(nr, rg)
    i2(nr, rg)
    i2x1(nr, rg)
    i2x2(nr, rg)
    i2x2d2(nr, rg)
    i2x2d1(nr, rg)
    i2x2div(nr, rg)
    i2x2laplh(nr, rg)
    i4(nr, rg)
    i4x4(nr, rg)
    i4x4laplh(nr, rg)
    i4x4lapl2h(nr, rg)
    qid(nr, rg)
