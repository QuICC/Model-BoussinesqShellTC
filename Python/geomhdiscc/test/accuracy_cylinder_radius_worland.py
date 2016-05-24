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
    lhs = transf.torspec(x_to_phys(res_expr,grid), m, op.shape[0])
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torspec(t, m, op.shape[0])
    err = np.abs(rhs[0:-(1+q)] - sol[q:-1])
    relerr = err/(1.0 + np.abs(sol[q:-1]))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err.T)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr.T)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, m, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    rhs = transf.torspec(x_to_phys(res_expr,grid), m, opA.shape[0])
    rhs = rhs[0:opA.shape[0]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = np.reshape(lhs, (lhs.shape[0],1))
    sol = transf.torspec(x_to_phys(sol_expr,grid), m, opA.shape[0])
    sol = sol[0:opA.shape[0]]
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err.T)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr.T)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

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
        print("\tTest for m = " + str(m))
        A = geo.i2(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,1,2)])
        ssol = 4**2*sy.integrate(sy.integrate(sphys*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 1)

def i2laplh(nr, rg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')

    print("\tForward:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,2*nr,2)])
        ssol = 4**2*sy.integrate(sy.integrate(sphys*x,x)*x,x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 1)

    print("\tbc = 10:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, {0:10}).tocsr()
        B = geo.i2(nr, m, {0:0}).tocsr()
        ssol = sy.expand((x**2-1)*np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)]))
        sphys = sy.expand(sy.diff(ssol,x,x) + sy.diff(ssol,x)/x - m**2*ssol/x**2)
        test_backward_tau(A, B, m, sphys, ssol, rg)

    print("\tbc = 11:")
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i2laplh(nr, m, {0:11}).tocsr()
        B = geo.i2(nr, m, {0:0}).tocsr()
        ssol = sy.expand((x**2-1)**2*np.sum([np.random.ranf()*x**(i+m) for i in np.arange(0,2*nr,2)]))
        sphys = sy.expand(sy.diff(ssol,x,x) + sy.diff(ssol,x)/x - m**2*ssol/x**2)
        test_backward_tau(A, B, m, sphys, ssol, rg)

def i4(nr, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    for i in range(0,2):
        m = np.random.randint(1, nr-1)
        m = m + (m+i)%2
        print("\tTest for m = " + str(m))
        A = geo.i4(nr, m, geo.radbc.no_bc())
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        ssol = 4**4*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sphys*x,x)*x,x)*x)*x)*x**m
        sphys = x**m*sphys
        test_forward(A, m, sphys, ssol, rg, 2)

def i4laplh(nr, rg):
    """Accuracy test for i4laplh operator"""

    print("i4laplh:")
    x = sy.Symbol('x')
    print("\tForward:")
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
        sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,2*nr,2)])
        ssol = 4**6*sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sy.integrate(sphys*x,x)*x,x)*x)*x)*x)*x)*x**m
        sphys = x**m*sphys
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

def test_fft(n, m):
    """Test Worland transform"""

    nr = np.ceil(3.0*n/2. + 3.0*m/4.0 + 1.0)
    # Create physical space test function
    t = wb.worland_poly(0, m, nr)
    for i in range(1, n):
        t += np.random.ranf()*wb.worland_poly(i, m, nr)
    t = np.matrix(t).T
    
    # Compute spectral expansion
    s = transf.torspec(t, m, n)

    # Project spectral expansion to physical space
    st = transf.torphys(s, m, nr)

    # Print error for transform loop
    print(np.max(np.abs(st.T-t.T)))

    import scipy.fftpack as fft
    import scipy.linalg as linalg
    cheb2worl  = []
    for i in range(0, n):
        cheb2worl.append(fft.dct(wb.worland_poly(i, m, nr))/(2*nr))
    cheb2worl = np.matrix(cheb2worl)
    f = fft.dct(t.T).T/(2*nr)
    mat = cheb2worl[:,m//2:m//2+n].T
    print(cheb2worl[:,16].T)
    print(cheb2worl[:,17].T)
    fs = linalg.solve_triangular(mat, f[m//2:m//2+n])

    # Print error for FFT transform loop
    print(np.max(np.abs(s.T-fs.T)))


if __name__ == "__main__":
    # Set test parameters
    n = 32
    nr = int(np.ceil(3.0*n/2.0 + 3.0*n/4.0 + 1))
    print("Grid: " + str((n, nr)))
    rg = wb.worland_grid(nr)

    # run tests
#    test_worland(nr, 110)
#    test_fft(nr, 32)
#    zblk(nr, rg)
#    i2(n, rg)
    i2laplh(nr, rg)
#    i4(n, rg)
#    i4laplh(nr, rg)
#    i4lapl2h(nr, rg)
#    i6(nr, rg)
#    i6laplh(nr, rg)
#    i6lapl2h(nr, rg)
#    i6lapl3h(nr, rg)
