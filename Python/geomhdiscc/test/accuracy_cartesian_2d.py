"""Check accuracy for cartesian 2D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin
import matplotlib.pylab as pl

import geomhdiscc.transform.cartesian as transf
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d


def vis_error(err, title):
    """Visualize the error"""

    if np.max(err) > 10*np.spacing(1):
        print(err)
        pl.imshow(np.log10(np.abs(err)))
        pl.title(title)
        pl.colorbar()
        pl.show()

def xz_to_phys(expr, grid_x, grid_z):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z), expr)
    vx, vz = np.meshgrid(grid_x, grid_z, indexing = 'ij')
    return func(vx, vz)

def test_forward(op, res_expr, sol_expr, grid_x, grid_z, qx, qz):
    """Perform a forward operation test"""

    print("\tForward test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    mesh = xz_to_phys(res_expr, grid_x, grid_z)
    lhs = transf.tocheb2d(mesh)
    lhs = lhs.reshape(nx*nz, order = 'F')
    rhs = op*lhs
    rhs = rhs.reshape(nx,nz, order='F')
    mesh = xz_to_phys(sol_expr,grid_x, grid_z)
    sol = transf.tocheb2d(mesh)
    err = np.abs(rhs - sol)
    vis_error(err, 'Forward error')
    print("\t\tMax forward error: " + str(np.max(err)))

def test_backward_tau(opA, opB, res_expr, sol_expr, grid_x, grid_z):
    """Perform a tau backward operation test"""

    print("\tBackward tau test")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    rhs = xz_to_phys(res_expr, grid_x, grid_z)
    rhs = transf.tocheb2d(rhs)
    rhs = rhs.reshape(nx*nz, order = 'F')
    lhs = spsplin.spsolve(opA,opB*rhs)
    lhs = lhs.reshape(nx,nz, order='F')
    sol = xz_to_phys(sol_expr, grid_x, grid_z)
    sol = transf.tocheb2d(sol)
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

def i2j2(nx,nz, xg, zg):
    """Accuracy test for i2j2 operator"""

    print("i2j2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

def i2j2d2d2(nx,nz, xg, zg):
    """Accuracy test for i2j2d2d2 operator"""

    print("i2j2d2d2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    A = c2d.i2j2d2d2(nx,nz, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sy.diff(sphys,x,x),z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 20, 20")
    k = np.random.ranf()*nx
    A = c2d.i2j2d2d2(nx,nz, {'x':{0:20}, 'z':{0:20}}).tocsr()
    B = c2d.i2j2(nx,nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(sy.diff(ssol,x,x),z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i2j0laplh(nx,nz, xg, zg):
    """Accuracy test for i2j0laplh operator"""

    print("i2j0laplh:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j0laplh(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys)
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 20, 0")
    k = np.random.ranf()*nx
    A = c2d.i2j0laplh(nx, nz, k, {'x':{0:20}, 'z':{0:0}}).tocsr()
    B = c2d.i2j0(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol)
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i2j2lapl(nx,nz, xg, zg):
    """Accuracy test for i2j2lapl operator"""

    print("i2j2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j2lapl(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    ssol = sy.integrate(ssol,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 20, 20")
    k = np.random.ranf()*nx
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

def i4j4lapl2(nx,nz, xg, zg):
    """Accuracy test for i4j4lapl2 operator"""

    print("i4j4lapl2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, c2d.c2dbc.no_bc())
    sphys = np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)]) for j in np.arange(0,nz,1)])
    ssol = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    ssol = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    ssol = sy.integrate(ssol,x,x,x,x)
    ssol = sy.expand(ssol)
    ssol = sy.integrate(ssol,z,z,z,z)
    ssol = sy.expand(ssol)
    test_forward(A, sphys, ssol, xg, zg, 2, 2)

    print("\tbc = 40, 40 (1)")
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}).tocsr()
    B = c2d.i4j4(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    sphys = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

    print("\tbc = 40, 40 (2)")
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}).tocsr()
    B = c2d.i4j4lapl(nx, nz, k, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    test_backward_tau(A, B, sphys, ssol, xg, zg)

if __name__ == "__main__":
    # Set test parameters
    nx = 32
    nz = 32
    xg = transf.grid(nx)
    zg = transf.grid(nz)

    # run hardcoded operator tests
    print('Hard coded exact operators')
    i2j2(nx, nz, xg, zg)
    i2j2d2d2(nx, nz, xg, zg)
    i2j0laplh(nx, nz, xg, zg)
    i2j2lapl(nx, nz, xg, zg)
    i4j4lapl2(nx, nz, xg, zg)
