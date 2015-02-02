"""Check accuracy for radial direction in a spherical shell"""

from __future__ import division
from __future__ import unicode_literals

import sympy as sy
import mpmath
import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl

import geomhdiscc.transform.shell as transf
import geomhdiscc.geometry.spherical.shell_radius as shell

import mpmath
mpmath.mp.dps = 200


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.lambdify(x, expr, 'mpmath')
    return func(grid)

def test_forward(op, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    x = sy.Symbol('x')
    lhs = transf.torcheb(x_to_phys(res_expr,grid))
    lhs = lhs[0:op.shape[0]]
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.torcheb(t)
    sol = sol[0:op.shape[0]]
    err = np.abs(rhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err[q:]) > 10*np.spacing(1):
        print(err)
    print("\t\tMax forward error: " + str(np.max(err[q:])))
    if np.max(relerr[q:]) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax forward relative error: " + str(np.max(relerr[q:])))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid))
    rhs = rhs[0:opA.shape[0]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torcheb(x_to_phys(sol_expr,grid))
    sol = sol[0:opA.shape[0]]
    err = np.abs(lhs - sol)
    relerr = err/(1.0 + np.abs(sol))
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))
    if np.max(relerr) > 10*np.spacing(1):
        print(relerr)
    print("\t\tMax tau backward relative error: " + str(np.max(relerr)))

def test_backward_galerkin(opA, opB, opS, res_expr, sol_expr, grid):
    """Perform a galerkin backward operation test"""

    print("\tBackward galkerin test")
    x = sy.Symbol('x')
    rhs = transf.torcheb(x_to_phys(res_expr,grid))
    rhs = rhs[0:opB.shape[1]]
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.torcheb(x_to_phys(sol_expr,grid))
    sol = sol[0:opB.shape[1]]
    tmpOp = opS[0:opS.shape[1],:]
    tmp = spsplin.spsolve(tmpOp, sol[0:opS.shape[1]])
    print(np.abs(opS*tmp - sol))
    err = np.abs(opS*lhs - sol)
    if np.max(err) > 10*np.spacing(1):
        print(err)
    print("\t\tMax tau backward error: " + str(np.max(err)))

def zblk(nr, a, b, rg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    A = shell.zblk(nr, 0, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = 0
    test_forward(A, sphys, ssol, rg, 0)

def d1(nr, a, b, rg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    A = shell.d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.diff(sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

def d2(nr, a, b, rg):
    """Accuracy test for d^2 operator"""

    print("d2:")
    x = sy.Symbol('x')
    A = shell.d2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.diff(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def divx(nr, a, b, rg):
    """Accuracy test for 1/x operator"""

    print("1/x:")
    print("\t Forward: (full spectrum):")
    A = shell.x1(nr, a, b, shell.radbc.no_bc()).tocsr()
    lhs = np.ones(rg.shape[0])
    sol = transf.torcheb(transf.torphys(lhs)/rg)
    rhs = spsplin.spsolve(A,lhs)
    print(sol-rhs)
    print("\t\tMax error for A^n: " + str(np.max(np.abs(sol-rhs))))
    pl.semilogy(abs(A*sol - np.ones(rg.shape[0]))+1e-16)
    pl.show()

def divx2(nr, a, b, rg):
    """Accuracy test for 1/x^2 operator"""

    print("1/x^2:")
    print("\t Forward: (full spectrum):")
    A = shell.x1(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = shell.x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    lhs = np.ones(rg.shape[0])
    sol = transf.torcheb(transf.torphys(lhs)/(rg*rg))
    tmp = spsplin.spsolve(A,lhs)
    rhsA = spsplin.spsolve(A,tmp)
    rhsB = spsplin.spsolve(B,lhs)
    print("\t\tMax error for A^n: " + str(np.max(np.abs(sol-rhsA))))
    print("\t\tMax error for B: " + str(np.max(np.abs(sol-rhsB))))
    pl.semilogy(abs(B*sol - np.ones(rg.shape[0]))+1e-16)
    pl.show()

def divx4(nr, a, b, rg):
    """Accuracy test for 1/x^4 operator"""

    print("1/x^4:")
    print("\t Forward: (full spectrum):")
    A = shell.x1(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = shell.x4(nr, a, b, shell.radbc.no_bc()).tocsr()
    lhs = np.ones(rg.shape[0])
    sol = transf.torcheb(transf.torphys(lhs)/(rg*rg*rg*rg))
    tmp = spsplin.spsolve(A,lhs)
    tmp = spsplin.spsolve(A,tmp)
    tmp = spsplin.spsolve(A,tmp)
    rhsA = spsplin.spsolve(A,tmp)
    rhsB = spsplin.spsolve(B,lhs)
    print("\t\tMax error for A^n: " + str(np.max(np.abs(sol-rhsA))))
    print("\t\tMax error for B: " + str(np.max(np.abs(sol-rhsB))))
    pl.semilogy(abs(B*sol - np.ones(rg.shape[0]))+1e-16)
    pl.show()

def x1(nr, a, b, rg):
    """Accuracy test for x1 operator"""

    print("x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full chebyshev):")
    A = shell.x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def x2(nr, a, b, rg):
    """Accuracy test for x^2 operator"""

    print("x2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.x2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full Chebyshev):")
    A = shell.x2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def x4(nr, a, b, rg):
    """Accuracy test for x^4 operator"""

    print("x4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.x4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*x*x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

    print("\t Forward: (full Chebyshev):")
    A = shell.x4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = (x*x*x*x*sphys)
    test_forward(A, sphys, ssol, rg, 0)

def i1(nr, a, b, rg):
    """Accuracy test for i1 operator"""

    print("i1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full Chebyshev):")
    A = shell.i1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.simplify(sy.integrate(sphys,x))
    test_forward(A, sphys, ssol, rg, 1)

def i1x1(nr, a, b, rg):
    """Accuracy test for i1x1 operator"""

    print("i1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(x*sphys,x)
    test_forward(A, sphys, ssol, rg, 1)

    print("\t Forward: (full Chebyshev):")
    A = shell.i1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.simplify(sy.integrate(x*sphys,x))
    test_forward(A, sphys, ssol, rg, 1)

def i2(nr, a, b, rg):
    """Accuracy test for i2 operator"""

    print("i2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2(nr, a, b, shell.radbc.no_bc()).tocsr()
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x1(nr, a, b, rg):
    """Accuracy test for i2x1 operator"""

    print("i2x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(x*sphys,x,x)
    test_forward(A, sphys, ssol, rg, 0)

def i2x2d1(nr, a, b, rg):
    """Accuracy test for i2x2d1 operator"""

    print("i2x2d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2x2d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x2d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x1d1x1(nr, a, b, rg):
    """Accuracy test for i2x1d1x1 operator"""

    print("i2x1d1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2x1d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x1d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2(nr, a, b, rg):
    """Accuracy test for i2x2 operator"""

    print("i2x2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2x2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**2)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x3(nr, a, b, rg):
    """Accuracy test for i2x3 operator"""

    print("i2x3:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i2x3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**3)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i2x2lapl(nr, a, b, rg):
    """Accuracy test for zblk operator"""

    print("i2x2lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x2lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**2*sy.diff(sphys,x,x) + 2*x*sy.diff(sphys,x) - l*(l+1)*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\tbc = 20:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:20}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand((x-(a+b))*(x-(-a+b))*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -20:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:-20, 'rt':2}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-20, 'rt':2}).tocsr()
    ssol = (x-(a+b))*(x-(-a+b))*np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr-2,1)])
#    ssol = sy.expand((x-(a+b))*(x-(-a+b))*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 21:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:21}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -21:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:-21, 'rt':2}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-21, 'rt':2}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 22:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:22, 'c':{'a':a, 'b':b}}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -22:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x2lapl(nr, l, a, b, {0:-22, 'rt':2, 'c':{'a':a, 'b':b}}).tocsr()
    B = shell.i2x2(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[2:,:]
    S = shell.radbc.stencil(nr, {0:-22, 'rt':2, 'c':{'a':a, 'b':b}}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) + 2*sy.diff(ssol,x)/x - l*(l+1)*ssol/x**2)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

def i2x3lapl(nr, a, b, rg):
    """Accuracy test for i2x3lapl operator"""

    print("i2x3lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i2x3lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x) + 2*x**2*sy.diff(sphys,x) - l*(l+1)*sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

    print("\t Forward: (full Chebyshev):")
    A = shell.i2x3lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x) + 2*x**2*sy.diff(sphys,x) - l*(l+1)*sphys*x)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, rg, 2)

def i4(nr, a, b, rg):
    """Accuracy test for i4 operator"""

    print("i4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4(nr, a, b, shell.radbc.no_bc()).tocsr()
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x1(nr, a, b, rg):
    """Accuracy test for i4x1 operator"""

    print("i4x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x1d1x1(nr, a, b, rg):
    """Accuracy test for i4x1d1x1 operator"""

    print("i4x1d1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x1d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x1d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x3(nr, a, b, rg):
    """Accuracy test for i4x3 operator"""

    print("i4x3:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x3(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4d1(nr, a, b, rg):
    """Accuracy test for i4x4d1 operator"""

    print("i4x4d1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x4d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x4d1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x3d1x1(nr, a, b, rg):
    """Accuracy test for i4x3d1x1 operator"""

    print("i4x3d1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x3d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x3d1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys*x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x3d2(nr, a, b, rg):
    """Accuracy test for i4x3d2 operator"""

    print("i4x3d2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x3d2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x3d2(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4(nr, a, b, rg):
    """Accuracy test for i4x4 operator"""

    print("i4x4:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x4(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(sphys*x**4)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4laplrd1x1(nr, a, b, rg):
    """Accuracy test for i4x4laplrd1x1 operator"""

    print("i4x4laplrd1x1:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    A = shell.i4x4laplrd1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x) + 3.0*x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x4laplrd1x1(nr, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x) + 3.0*x**3*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl(nr, a, b, rg):
    """Accuracy test for i4x4lapl operator"""

    print("i4x4lapl:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x4lapl(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x) + 2*x**3*sy.diff(sphys,x) - l*(l+1)*x**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

def i4x4lapl2(nr, a, b, rg):
    """Accuracy test for i4x4lapl2 operator"""

    print("i4x4lapl2:")
    print("\t Forward: (polynomial):")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl2(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\t Forward: (full Chebyshev):")
    A = shell.i4x4lapl2(nr, l, a, b, shell.radbc.no_bc())
    sphys = np.sum([(-1.0)**np.round(np.random.ranf())*sy.chebyshevt(int(i),(x-b)/a) for i in np.arange(0,nr,1)])
    ssol = sy.expand(x**4*sy.diff(sphys,x,x,x,x) + 4*x**3*sy.diff(sphys,x,x,x) - 2*l*(l+1)*x**2*sy.diff(sphys,x,x) + (l-1)*l*(l+1)*(l+2)*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, rg, 4)

    print("\tbc = 40:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl2(nr, l, a, b, {0:40}).tocsr()
    B = shell.i4x4(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -40:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl2(nr, l, a, b, {0:-40, 'rt':4}).tocsr()
    B = shell.i4x4(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[4:,:]
    S = shell.radbc.stencil(nr, {0:-40, 'rt':4}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**2*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

    print("\tbc = 41:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl2(nr, l, a, b, {0:41}).tocsr()
    B = shell.i4x4(nr, a, b, shell.radbc.no_bc()).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**3*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-6,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_tau(A, B, sphys, ssol, rg)

    print("\tbc = -41:")
    x = sy.Symbol('x')
    l = np.random.randint(1, nr)
    A = shell.i4x4lapl2(nr, l, a, b, {0:-41, 'rt':4}).tocsr()
    B = shell.i4x4(nr, a, b, shell.radbc.no_bc()).tocsr()
    B = B[4:,:]
    S = shell.radbc.stencil(nr, {0:-41, 'rt':4}).tocsr()
    ssol = sy.expand(((x-(a+b))*(x-(-a+b)))**3*np.sum([np.random.ranf()*x**(i) for i in np.arange(0,nr-3,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + 4*sy.diff(ssol,x,x,x)/x - 2*l*(l+1)*sy.diff(ssol,x,x)/x**2 + (l-1)*l*(l+1)*(l+2)*ssol/x**4)
    test_backward_galerkin(A, B, S, sphys, ssol, rg)

def qid(nr, a, b, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    A = shell.qid(nr, 3, shell.radbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nr,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)

if __name__ == "__main__":
    # Set test parameters
    nr = 14
    a, b = shell.linear_r2x(1.0, 0.35)
    rg = transf.rgrid(2*nr, a, b)

    # run tests
    #zblk(nr, a, b, rg)
#    d1(nr, a, b, rg)
#    d2(nr, a, b, rg)
#    x1(nr, a, b, rg)
#    x2(nr, a, b, rg)
#    x4(nr, a, b, rg)
#    i1(nr, a, b, rg)
#    i1x1(nr, a, b, rg)
#    i2(nr, a, b, rg)
#    i2x1(nr, a, b, rg)
#    i2x1d1x1(nr, a, b, rg)
#    i2x2d1(nr, a, b, rg)
#    i2x2(nr, a, b, rg)
#    i2x3(nr, a, b, rg)
    i2x2lapl(nr, a, b, rg)
#    i2x3lapl(nr, a, b, rg)
#    i4(nr, a, b, rg)
#    i4x1(nr, a, b, rg)
#    i4x1d1x1(nr, a, b, rg)
#    i4x3d2(nr, a, b, rg)
#    i4x3(nr, a, b, rg)
#    i4x3d1x1(nr, a, b, rg)
#    i4x4d1(nr, a, b, rg)
#    i4x4(nr, a, b, rg)
#    i4x4laplrd1x1(nr, a, b, rg)
#    i4x4lapl(nr, a, b, rg)
    i4x4lapl2(nr, a, b, rg)
#    qid(nr, a, b, rg)

#    divx(nr, a, b, rg)
#    divx2(nr, a, b, rg)
#    divx4(nr, a, b, rg)
