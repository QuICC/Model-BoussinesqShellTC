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
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.zblk(nr, 0, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = 0.0
        test_forward(A, l%2, sphys, ssol, rg, 0)

def d1(nr, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.diff(sphys,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)
def i2x1(nr, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2x1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)

def i2x2d1(nr, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2x2d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 1)


def i2x2(nr, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2x2(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**2)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)


def i2x2lapl(nr, rg):
    """Accuracy test for zblk operator"""

    print("i2x2lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i2x2lapl(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 1)

def i4x3(nr, rg):
    """Accuracy test for i4x3 operator"""

    print("i4x3:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4x3(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**3*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)

def i4x4d1(nr, rg):
    """Accuracy test for i4x4d1 operator"""

    print("i4x4d1:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4x4d1(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x))
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, (l%2,(l+1)%2), sphys, ssol, rg, 2)

def i4x4(nr, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4x4(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def i4x4lapl(nr, rg):
    """Accuracy test for i4x4lapl operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4x4lapl(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def i4x4lapl2(nr, rg):
    """Accuracy test for i4x4lapl2 operator"""

    print("i4x4lapl2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        print("\tTest for l = " + str(l))
        A = sphere.i4x4lapl2(nr, l, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(l%2,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, l%2, sphys, ssol, rg, 2)

def qid(nr, rg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    for i in range(0,2):
        l = np.random.randint(1, nr-1)
        l = l + (l+i)%2
        A = sphere.qid(nr, l, 3, sphere.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(l%2,2*nr,2)])
        ssol = sphys
        test_forward(A, l%2, sphys, ssol, rg, 3)


if __name__ == "__main__":
    # Set test parameters
    nr = 20
    rg = transf.rgrid(nr)

    # run tests
    #zblk(nr, rg)
    d1(nr, rg)
    i2x1(nr, rg)
    i2x2d1(nr, rg)
    i2x2(nr, rg)
    i2x2lapl(nr, rg)
    i4x3(nr, rg)
    i4x4d1(nr, rg)
    i4x4(nr, rg)
    i4x4lapl(nr, rg)
    i4x4lapl2(nr, rg)
    qid(nr, rg)
