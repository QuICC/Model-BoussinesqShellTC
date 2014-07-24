"""Check accuracy for cartesian 1D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl

import geomhdiscc.transform.cartesian as transf
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_generic_1d as cg1d


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        pl.semilogy(np.abs(err))
        pl.title(title)
        pl.show()

def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

def test_forward(op, res_expr, sol_expr, grid, q):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    lhs = transf.tocheb(x_to_phys(res_expr,grid))
    rhs = op*lhs
    t = x_to_phys(sol_expr,grid)
    sol = transf.tocheb(t)
    pl.semilogy(np.abs(sol))
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err[q:])))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid):
    """Perform a tau backward operation test"""

    print("\tBackward tau test")
    x = sy.Symbol('x')
    rhs = transf.tocheb(x_to_phys(res_expr,grid))
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.tocheb(x_to_phys(sol_expr,grid))
    err = np.abs(lhs - sol)
    vis_error(err, 'Tau backward error')
    print("\t\tMax tau backward error: " + str(np.max(err)))

def test_backward_galerkin(opA, opB, opS, res_expr, sol_expr, grid):
    """Perform a galerkin backward operation test"""

    print("\tBackward galkerin test")
    x = sy.Symbol('x')
    rhs = transf.tocheb(x_to_phys(res_expr,grid))
    lhs = spsplin.spsolve(opA,opB*rhs)
    sol = transf.tocheb(x_to_phys(sol_expr,grid))
    err = np.abs(opS*lhs - sol)
    vis_error(err, 'Galerkin backward error')
    print("\t\tMax galerkin backward error: " + str(np.max(err)))

def zblk(nx, xg):
    """Accuracy test for zblk operator"""

    print("zblk:")
    x = sy.Symbol('x')
    A = c1d.zblk(nx)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = 0
    test_forward(A, sphys, ssol, xg, 0)

def d1(nx, xg):
    """Accuracy test for d1 operator"""

    print("d1:")
    x = sy.Symbol('x')
    A = c1d.d1(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.diff(sphys,x)
    test_forward(A, sphys, ssol, xg, 1)

def d2(nx, xg):
    """Accuracy test for d2 operator"""

    print("d2:")
    x = sy.Symbol('x')
    A = c1d.d2(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.diff(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 2)

def d4(nx, xg):
    """Accuracy test for d4 operator"""

    print("d4:")
    x = sy.Symbol('x')
    A = c1d.d4(nx, c1d.c1dbc.no_bc())
    print(A.todense())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.diff(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

def laplh(nx, xg):
    """Accuracy test for laplh operator"""

    print("laplh:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.laplh(nx, k, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.diff(sphys,x,x) - k**2*sphys
    test_forward(A, sphys, ssol, xg, 2)

def lapl2h(nx, xg):
    """Accuracy test for lapl2h operator"""

    print("lapl2h:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.lapl2h(nx, k, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) - 2.0*k**2*sy.diff(sphys,x,x) + k**4*sphys)
    test_forward(A, sphys, ssol, xg, 4)

def i1(nx, xg):
    """Accuracy test for i1 operator"""

    print("i1:")
    x = sy.Symbol('x')
    A = c1d.i1(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sy.integrate(sphys,x)
    test_forward(A, sphys, ssol, xg, 1)

def i2(nx, xg):
    """Accuracy test for i2 operator"""

    print("i2:")
    x = sy.Symbol('x')
    A = c1d.i2(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.integrate(sphys,x,x)
    test_forward(A, sphys, ssol, xg, 2)

def i2d2(nx, xg):
    """Accuracy test for i2d2 operator"""

    print("i2d2:")
    x = sy.Symbol('x')
    A = c1d.i2d2(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)

    print("\tbc = 20")
    A = c1d.i2d2(nx, {0:20}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.diff(ssol,x,x)
    test_backward_tau(A, B, sphys, ssol, xg)

    print("\tbc = -20")
    A = c1d.i2d2(nx, {0:-20, 'r':2}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tolil()
    B = B[2:,:]
    S = c1d.c1dbc.stencil(nx, {0:-20})
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.diff(ssol,x,x)
    #test_backward_galerkin(A, B, S, sphys, ssol, xg)

def i2d1(nx, xg):
    """Accuracy test for i2d1 operator"""

    print("i2d1:")
    x = sy.Symbol('x')
    A = c1d.i2d1(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.integrate(sy.diff(sphys,x),x,x)
    test_forward(A, sphys, ssol, xg, 2)

def i2lapl(nx, xg):
    """Accuracy test for i2lapl operator"""

    print("i2lapl:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i2lapl(nx, k, l, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys - l**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)

    print("\tbc = 20")
    A = c1d.i2lapl(nx, k, l, {0:20}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

    print("\tbc = -20")
    A = c1d.i2lapl(nx, k, l, {0:-20, 'r':2}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    B = B[2:,:]
    S = c1d.c1dbc.stencil(nx, {0:-20})
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    #test_backward_galerkin(A, B, S, sphys, ssol, xg)

def i2laplh(nx, xg):
    """Accuracy test for i2laplh operator"""

    print("i2laplh:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i2laplh(nx, k, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    test_forward(A, sphys, ssol, xg, 2)

    print("\tbc = 20")
    A = c1d.i2laplh(nx, k, {0:20}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

    print("\tbc = -20")
    A = c1d.i2laplh(nx, k, {0:-20, 'r':2}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    B = B[2:,:]
    S = c1d.c1dbc.stencil(nx, {0:-20})
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    #test_backward_galerkin(A, B, S, sphys, ssol, xg)

def i4(nx, xg):
    """Accuracy test for i4 operator"""

    print("i4:")
    x = sy.Symbol('x')
    A = c1d.i4(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)])
    ssol = sy.integrate(sphys,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

def i4d2(nx, xg):
    """Accuracy test for i4d2 operator"""

    print("i4d2:")
    x = sy.Symbol('x')
    A = c1d.i4d2(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

def i4d4(nx, xg):
    """Accuracy test for i4d4 operator"""

    print("i4d4:")
    x = sy.Symbol('x')
    A = c1d.i4d4(nx, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

def i4lapl(nx, xg):
    """Accuracy test for i4lapl operator"""

    print("i4lapl:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl(nx, k, l, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys - l**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

    print("\tbc = 40")
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl(nx, k, l, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

def i4laplh(nx, xg):
    """Accuracy test for i4laplh operator"""

    print("i4laplh:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i4laplh(nx, k, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

    print("\tbc = 40")
    k = np.random.ranf()*nx
    A = c1d.i4laplh(nx, k, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

def i4lapl2(nx, xg):
    """Accuracy test for i4lapl2 operator"""

    print("i4lapl2:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl2(nx, k, l, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) + k**4*sphys + l**4*sphys - 2*k**2*sy.diff(sphys,x,x) - 2*l**2*sy.diff(sphys,x,x) + 2*k**2*l**2*sphys)
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

    print("\tbc = 40 (1)")
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl2(nx, k, l, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + k**4*ssol + l**4*ssol - 2*k**2*sy.diff(ssol,x,x) - 2*l**2*sy.diff(ssol,x,x) + 2*k**2*l**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

    print("\tbc = 40 (2)")
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl2(nx, k, l, {0:40}).tocsr()
    B = c1d.i4lapl(nx, k, l, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg)

    print("\tbc = 40 split (2)")
    A11 = c1d.i2lapl(nx, k, l, {0:20}).tolil()
    A12 = -c1d.i2(nx, c1d.c1dbc.no_bc()).tolil()
    A21 = c1d.zblk(nx, {0:21}).tolil()
    A22 = c1d.i2lapl(nx, k, l, {0:0}).tolil()
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B11 = c1d.zblk(nx, c1d.c1dbc.no_bc()).tolil()
    B12 = c1d.zblk(nx, c1d.c1dbc.no_bc()).tolil()
    B21 = c1d.zblk(nx, c1d.c1dbc.no_bc()).tolil()
    B22 = c1d.i2lapl(nx, k, l, c1d.c1dbc.no_bc()).tolil()
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    rhs = transf.tocheb(x_to_phys(sphys,xg))
    rhs = np.append(0*rhs, rhs)
    lhs = spsplin.spsolve(A,B*rhs)
    sol = transf.tocheb(x_to_phys(ssol,xg))
    sol = np.append(sol, 0*sol)
    err = np.abs(lhs - sol)
    if np.max(err[0:nx]) > 10*np.spacing(1):
        print(err[0:nx])
    print("\t\tMax forward error: " + str(np.max(err[0:nx])))

    print("\tbc = 40 split (trunc) (2)")
    A12 = A12[:,0:-2]
    A22 = A22[:,0:-2]
    A21 = A21[0:-2,:]
    A22 = A22[0:-2,:]
    A = spsp.bmat([[A11,A12],[A21,A22]]).tocsr()
    B12 = B12[:,0:-2]
    B22 = B22[:,0:-2]
    B21 = B21[0:-2,:]
    B22 = B22[0:-2,:]
    B = spsp.bmat([[B11,B12],[B21,B22]]).tocsr()
    rhs = transf.tocheb(x_to_phys(sphys,xg))
    rhs = np.append(0*rhs, rhs[0:-2])
    lhs = spsplin.spsolve(A,B*rhs)
    sol = transf.tocheb(x_to_phys(ssol,xg))
    sol = np.append(sol, 0*sol[0:-2])
    err = np.abs((lhs - sol))
    if np.max(err[0:nx]) > 10*np.spacing(1):
        print(err[0:nx])
    print("\t\tMax forward error: " + str(np.max(err[0:nx])))

def i4lapl2h(nx, xg):
    """Accuracy test for i4lapl2h operator"""

    print("i4lapl2h:")
    x = sy.Symbol('x')
    k = np.random.ranf()*nx
    A = c1d.i4lapl2h(nx, k, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sy.expand(sy.diff(sphys,x,x,x,x) + k**4*sphys - 2*k**2*sy.diff(sphys,x,x))
    ssol = sy.integrate(ssol,x,x,x,x)
    test_forward(A, sphys, ssol, xg, 4)

    print("\tbc = 40")
    k = np.random.ranf()*nx
    A = c1d.i4lapl2h(nx, k, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + k**4*ssol - 2*k**2*sy.diff(ssol,x,x))
    test_backward_tau(A, B, sphys, ssol, xg)

def qid(nx, xg):
    """Accuracy test for qid operator"""

    print("qid:")
    x = sy.Symbol('x')
    A = c1d.qid(nx, 3, c1d.c1dbc.no_bc())
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
    ssol = sphys
    test_forward(A, sphys, ssol, xg, 3)

def genxp(nx, q, p, xg, ntrunc = -1):
    """Accuracy test for genxp operator"""

    print("genx"+ str(p) + ":")
    x = sy.Symbol('x')
    expr = 'x**' + str(p)
    A = cg1d.generic(nx, q, expr, x, c1d.c1dbc.no_bc(), 1.0, ntrunc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sphys*x**p
    test_forward(A, sphys, ssol, xg, 1)

def genlinxp(nx, q, p, xg, ntrunc = -1):
    """Accuracy test for gen_linxp operator"""

    print("genlinx"+ str(p) + ":")
    x = sy.Symbol('x')
    expr = '(1+x)**' + str(p)
    print(expr)
    A = cg1d.generic(nx, q, expr, x, c1d.c1dbc.no_bc(), 1.0, ntrunc)
    sphys = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
    ssol = sphys*(1+x)**p
    test_forward(A, sphys, ssol, xg, 1)


if __name__ == "__main__":
    # Set test parameters
    nx = 30
    xg = transf.grid(nx)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    #zblk(nx, xg)
    d1(nx, xg)
    d2(nx, xg)
    d4(nx, xg)
    laplh(nx, xg)
    lapl2h(nx, xg)
    i1(nx, xg)
    i2(nx, xg)
    i2d1(nx, xg)
    i2d2(nx, xg)
    i2lapl(nx, xg)
    i2laplh(nx, xg)
    i4(nx, xg)
    i4d2(nx, xg)
    i4d4(nx, xg)
    i4lapl(nx, xg)
    i4laplh(nx, xg)
    i4lapl2(nx, xg)
    i4lapl2h(nx, xg)
    qid(nx, xg)

    # run generic operator tests
    print('Generic operator generator')
    genxp(nx, 0, 1, xg)
    genxp(nx, 0, 2, xg)
    genxp(nx, 0, 3, xg)
    genxp(nx, 0, 4, xg)
    genxp(nx, 0, 9, xg)
    genlinxp(nx, 0, 1, xg)
    genlinxp(nx, 0, 2, xg)
    genlinxp(nx, 0, 1.3, xg)
    genlinxp(nx, 0, 1.3, xg, 10)
