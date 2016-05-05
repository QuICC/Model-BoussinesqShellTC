"""Check accuracy for radial direction in a cylinder"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import scipy.special as special

import geomhdiscc.transform.cylinder_worland as transf
import geomhdiscc.geometry.cylindrical.cylinder_radius_worland as geo
import geomhdiscc.geometry.worland.worland_basis as wb


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

def test_forward(op, m, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.torspec(x_to_phys(res_expr,grid), m, nr)
    print("LHS:")
    print(lhs.T)
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torspec(t, m, nr)
    print("solution:")
    print(sol.T)
    print("computation:")
    print(rhs.T)
    print("ratio:")
    print(rhs.T/sol.T)
    print(rhs[0:-2].T/sol[1:-1].T)
    print("-----------------------------------------")
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def zblk(nr, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.zblk(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = 0
        test_forward(A, m, sphys, ssol, rg, 0)

def i2(nr, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        m = 3
        print("\tTest for m = " + str(m))
        A = geo.i2(nr, m, geo.radbc.no_bc())
        #sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,1,2)])
        sphys = x**m*wb.worland_norm(0,m)
        ssol = 16*sy.integrate(sy.integrate(sphys*x**(1-m),x)*x,x)*x**m
        test_forward(A, m, sphys, ssol, rg, 1)

def i2laplh(nr, rg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(sphys*x)
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, m, sphys, ssol, rg, 1)

def i4(nr, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x*sphys)
        ssol = sy.expand(x*sy.diff(ssol,x))
        ssol = sy.integrate(ssol,x,x)
        test_forward(A, m, sphys, ssol, rg, 2)

def i4laplh(nr, rg):
    """Accuracy test for i4laplh operator"""

    print("i4laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(sphys*x**4)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 2)

def i4lapl2h(nr, rg):
    """Accuracy test for i4lapl2h operator"""

    print("i4x4lapl:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4lapl2h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x) + x**3*sy.diff(sphys,x) - m**2*x**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 2)

def i6(nr, rg):
    """Accuracy test for i6 operator"""

    print("i6:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 3)

def i6laplh(nr, rg):
    """Accuracy test for i6laplh operator"""

    print("i6laplh:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 3)

def i6lapl2h(nr, rg):
    """Accuracy test for i6lapl2h operator"""

    print("i6lapl2h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6lapl2h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 3)

def i6lapl3h(nr, rg):
    """Accuracy test for i6lapl3h operator"""

    print("i6lapl3h:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i6lapl3h(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)])
        ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 2*x**3*sy.diff(sphys,x,x,x) - (1+2*m**2)*x**2*sy.diff(sphys,x,x) + (1+2*m**2)*x*sy.diff(sphys,x) + (m**2 - 4)*m**2*sphys)
        ssol = sy.integrate(ssol,x,x,x,x)
        test_forward(A, m, sphys, ssol, rg, 3)

def test_worland(nr, m):
    """Test Worland transform"""

    # Create physical space test function
    t = wb.worland_poly(0, m, nr)
    for i in range(1, nr-m):
        t += np.random.ranf()*wb.worland_poly(i, m, nr)
    t = np.matrix(t).T
    
    # Compute spectral expansion
    s = transf.torspec(t, m, nr-m)

    # Project spectral expansion to physical space
    st = transf.torphys(s, m, nr)

    # Print error for transform loop
    print(np.abs(st.T-t.T))


if __name__ == "__main__":
    # Set test parameters
    nr = 16
    rg = wb.worland_grid(np.ceil(5*nr))

    # run tests
#    test_worland(nr, 110)
#    zblk(nr, rg)
    i2(nr, rg)
#    i2laplh(nr, rg)
#    i4(nr, rg)
#    i4laplh(nr, rg)
#    i4lapl2h(nr, rg)
#    i6(nr, rg)
#    i6laplh(nr, rg)
#    i6lapl2h(nr, rg)
#    i6lapl3h(nr, rg)
