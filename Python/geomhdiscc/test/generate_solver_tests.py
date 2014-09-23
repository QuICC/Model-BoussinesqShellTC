"""Generate test problem matrices for sparse solver tests"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy as sy
import scipy.io as io
if True:
    import matplotlib.pylab as pl
    has_error_plot = True
else:
    has_error_plot = False

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.geometry.cartesian.cartesian_3d as c3d
import geomhdiscc.transform.cartesian as transf


def x_to_phys(expr, grid):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    func = sy.utilities.lambdify(x, expr)
    return func(grid)

def xz_to_phys(expr, grid_x, grid_z):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z), expr)
    vx, vz = np.meshgrid(grid_x, grid_z, indexing = 'ij')
    return func(vx, vz)

def xyz_to_phys(expr, grid_x, grid_z, grid_y):
    """Convert sympy expression to grid values"""

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    func = sy.utilities.lambdify((x, z, y), expr)
    vx, vz, vy = np.meshgrid(grid_x, grid_z, grid_y, indexing = 'ij')
    return func(vx, vz, vy)

def write1DTest(opA, opB, res_expr, sol_expr, grid, name):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    nx = len(grid)
    rhs = transf.tocheb(x_to_phys(res_expr,grid))
    sol = transf.tocheb(x_to_phys(sol_expr,grid))
    io.mmwrite(name + '_A_' +  str(nx) + '.mtx', opA)
    io.mmwrite(name + '_rhs_' +  str(nx) + '.mtx', (opB*rhs).reshape(nx, 1))
    io.mmwrite(name + '_sol_' +  str(nx) + '.mtx', sol.reshape(nx, 1))

def write2DTest(opA, opB, res_expr, sol_expr, grid_x, grid_z, name):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    z = sy.Symbol('z')
    nx = len(grid_x)
    nz = len(grid_x)
    rhs = xz_to_phys(res_expr, grid_x, grid_z)
    rhs = transf.tocheb2d(rhs)
    rhs = rhs.reshape(nx*nz, order = 'F')
    sol = xz_to_phys(sol_expr, grid_x, grid_z)
    sol = transf.tocheb2d(sol)
    io.mmwrite(name + '_A_' +  str(nx*nz) + '.mtx', opA)
    io.mmwrite(name + '_rhs_' + str(nx*nz) + '.mtx', (opB*rhs).reshape(nx*nz, 1, order = 'F'))
    io.mmwrite(name + '_sol_' + str(nx*nz) + '.mtx', sol.reshape(nx*nz, 1, order = 'F'))

def write3DTest(opA, opB, res_expr, sol_expr, grid_x, grid_y, grid_z, name):
    """Perform a tau backward operation test"""

    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')
    nx = len(grid_x)
    ny = len(grid_y)
    nz = len(grid_x)
    rhs = xyz_to_phys(res_expr, grid_x, grid_z, grid_y)
    rhs = transf.tocheb3d(rhs)
    rhs = rhs.reshape(nx*ny*nz, order='F')
    sol = xyz_to_phys(sol_expr, grid_x, grid_z, grid_y)
    sol = transf.tocheb3d(sol)
    io.mmwrite(name + '_A_' + str(nx*nz*nz) + '.mtx', opA)
    io.mmwrite(name + '_rhs_' + str(nx*nz*nz) + '.mtx', (opB*rhs).reshape(nx*ny*nz,1,order = 'F'))
    io.mmwrite(name + '_sol_' + str(nx*nz*nz) + '.mtx', sol.reshape(nx*ny*nz,1,order = 'F'))

def i2lapl(nx, xg):
    """Accuracy test for i2lapl operator"""

    print("i2lapl:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i2lapl(nx, k, l, {0:20}).tocsr()
    B = c1d.i2(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol - l**2*ssol)
    write1DTest(A, B, sphys, ssol, xg, 'laplacian1D')

def i4lapl2(nx, xg):
    """Accuracy test for i4lapl2 operator"""

    print("i4lapl2:")
    x = sy.Symbol('x')
    k, l = np.random.rand(2)*nx
    A = c1d.i4lapl2(nx, k, l, {0:40}).tocsr()
    B = c1d.i4(nx, c1d.c1dbc.no_bc()).tocsr()
    ssol = sy.expand((1.0 - x**2)**2*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]))
    sphys = sy.expand(sy.diff(ssol,x,x,x,x) + k**4*ssol + l**4*ssol - 2*k**2*sy.diff(ssol,x,x) - 2*l**2*sy.diff(ssol,x,x) + 2*k**2*l**2*ssol)
    write1DTest(A, B, sphys, ssol, xg, 'bilaplacian1D')

def i2j2lapl(nx,nz, xg, zg):
    """Accuracy test for i2j2lapl operator"""

    print("i2j2lapl:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i2j2lapl(nx, nz, k, {'x':{0:20}, 'z':{0:20}}).tocsr()
    B = c2d.i2j2(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    write2DTest(A, B, sphys, ssol, xg, zg, 'laplacian2D')

def i4j4lapl2(nx,nz, xg, zg):
    """Accuracy test for i4j4lapl2 operator"""

    print("i4j4lapl2:")
    x = sy.Symbol('x')
    z = sy.Symbol('z')
    k = np.random.ranf()*nx
    A = c2d.i4j4lapl2(nx, nz, k, {'x':{0:40}, 'z':{0:40}}).tocsr()
    B = c2d.i4j4(nx, nz, c2d.c2dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**2*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)]) for j in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.diff(ssol,x,x) - k**2*ssol + sy.diff(ssol,z,z))
    sphys = sy.expand(sy.diff(sphys,x,x) - k**2*sphys + sy.diff(sphys,z,z))
    write2DTest(A, B, sphys, ssol, xg, zg, 'bilaplacian2D')

def i2j2k2lapl(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i2j2k2lapl operator"""

    print("i2j2k2lapl:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')

    A = c3d.i2j2k2lapl(nx, ny, nz, {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'xy'}).tocsr()
    B = c3d.i2j2k2(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)*(1.0 - y**2)*(1.0 - z**2)*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)]) for j in np.arange(0,ny-2,1)]) for k in np.arange(0,nz-2,1)])
    sphys = sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z))
    write3DTest(A, B, sphys, ssol, xg, yg, zg, 'laplacian3D')

def i4j4k4lapl2(nx, ny, nz, xg, yg, zg):
    """Accuracy test for i4j4k4lapl2 operator"""

    print("i4j4k4lapl2:")
    x = sy.Symbol('x')
    y = sy.Symbol('y')
    z = sy.Symbol('z')

    A = c3d.i4j4k4lapl2(nx, ny, nz, {'x':{0:41}, 'y':{0:41}, 'z':{0:40}, 'priority':'zsx'}).tocsr()
    B = c3d.i4j4k4(nx, ny, nz, c3d.c3dbc.no_bc()).tocsr()
    ssol = (1.0 - x**2)**3*(1.0 - y**2)**3*(1.0 - z**2)**2*np.sum([np.random.ranf()*z**k*np.sum([np.random.ranf()*y**j*np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-6,1)]) for j in np.arange(0,ny-6,1)]) for k in np.arange(0,nz-4,1)])
    sphys = sy.expand(sy.expand(sy.diff(ssol,x,x)) + sy.expand(sy.diff(ssol,y,y)) + sy.expand(sy.diff(ssol,z,z)))
    sphys = sy.expand(sy.expand(sy.diff(sphys,x,x)) + sy.expand(sy.diff(sphys,y,y)) + sy.expand(sy.diff(sphys,z,z)))
    write3DTest(A, B, sphys, ssol, xg, yg, zg, 'bilaplacian3D')

if __name__ == "__main__":
    # Set test parameters
    ns = [16, 32, 64, 128, 256, 512]
    for nn in ns:
        xg = transf.grid(nn)
        i2lapl(nn, xg)
        i4lapl2(nn, xg)

    ns = [8, 16, 24, 32, 48]
    for nn in ns:
        xg = transf.grid(nn)
        zg = transf.grid(nn)
        i2j2lapl(nn, nn, xg, zg)
        i4j4lapl2(nn, nn, xg, zg)

    ns = [8, 10, 12, 14, 16]
    for nn in ns:
        xg = transf.grid(nn)
        yg = transf.grid(nn)
        zg = transf.grid(nn)
        i2j2k2lapl(nn, nn, nn, xg, yg, zg)
        i4j4k4lapl2(nn, nn, nn, xg, yg, zg)
