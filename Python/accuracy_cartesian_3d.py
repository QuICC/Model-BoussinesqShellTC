"""Check accuracy of operations for cartesian 3D operators"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import sympy
import chebyshev_tools as ct
import cartesian_3d as c3d
import scipy.sparse as spsp

x = sympy.Symbol('x')
z = sympy.Symbol('z')

def x_to_phys(expr, grid):
   func = sympy.utilities.lambdify(x, expr)
   return func(grid)


def test_forward(op, res_expr, sol_expr, grid, q):
   print("NOTHING YET")
#   lhs = ct.tocheb(x_to_phys(res_expr,grid))
#   rhs = op*lhs
#   t = x_to_phys(sol_expr,grid)
#   sol = ct.tocheb(t)
#   err = np.abs(rhs - sol)
#   print(sol)
#   print(rhs)
#   print(err)
#   print(np.max(err[q:]))


nx = 20
ny = 20
nz = 20
xg = ct.grid(nx)
yg = ct.grid(ny)
zg = ct.grid(nz)

#
# Accuracy tests: zblk
#
print("zblk:")
A = c3d.zblk(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz,1)])
#ssolx = 0
#ssolz = 0
#test_forward(A, sphys, ssol, xg, 0)

#
# Accuracy tests: i2j2d2d2
#
print("i2j2k2d2d2d2:")
A = c3d.i2j2k2d2d2d2(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz,1)])
#ssolx = sympy.integrate(sympy.diff(sphysx,x,x),x,x)
#ssolz = sympy.integrate(sympy.diff(sphysz,z,z),z,z)
#test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i2j2
#
print("i2j2k2:")
A = c3d.i2j2k2(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-1,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz-1,1)])
#ssolx = sympy.integrate(sphysx,x,x)
#ssolz = sympy.integrate(sphysz,z,z)
#test_forward(A, sphys, ssol, xg, 1)

#
# Accuracy tests: i2j2lapl
#
print("i2j2k2lapl:")
k = np.random.ranf()*nx
A = c3d.i2j2k2lapl(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-2,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz-2,1)])
#ssol = sympy.integrate((sympy.diff(sphysx,x,x)*sphysz).expand(),x,x,z,z) + sympy.integrate((sympy.diff(sphysz,z,z)*sphysx).expand(),x,x,z,z) - k**2*sympy.integrate((sphysx*sphysz).expand(),x,x,z,z)
#test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i4j4
#
print("i4j4k4:")
A = c3d.i4j4k4(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz,1)])
#ssolx = sympy.integrate(sphysx,x,x,x,x)
#ssolz = sympy.integrate(sphysz,z,z,z,z)
#test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i4j4lapl
#
print("i4j4k4lapl:")
k = np.random.ranf()*nx
A = c3d.i4j4k4lapl(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz,1)])
#ssol = sympy.simplify((sympy.diff(sphysx,x,x)*sphysz) + (sympy.diff(sphysz,z,z)*sphysx) - k**2*(sphysx*sphysz))
#ssol = sympy.integrate(ssol.expand() ,x,x,x,x,z,z,z,z)
#ssol = sympy.integrate((sympy.diff(sphysx,x,x)*sphysz).expand() ,x,x,x,x,z,z,z,z) + sympy.integrate((sympy.diff(sphysz,z,z)*sphysx).expand() ,x,x,x,x,z,z,z,z) - k**2*sympy.integrate((sphysx*sphysz).expand() ,x,x,x,x,z,z,z,z)
#test_forward(A, sphys, ssol, xg, 2)

#
# Accuracy tests: i4j4lapl2
#
print("i4j4k4lapl2:")
k = np.random.ranf()*nx
A = c3d.i4j4k4lapl2(nx,ny,nz)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx-4,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz-4,1)])
#ssol = sympy.integrate(sympy.diff(sphysx,x,x,x,x)*sphysz ,x,x,x,x,z,z,z,z) + sympy.integrate(sympy.diff(sphysz,z,z,z,z)*sphysx ,x,x,x,x,z,z,z,z) + k**4*sympy.integrate(sphysx*sphysz ,x,x,x,x,z,z,z,z) - k**2*sympy.integrate(sympy.diff(sphysx, x,x)*sphysz ,x,x,x,x,z,z,z,z) - k**2*sympy.integrate(sympy.diff(sphysz, z,z)*sphysx ,x,x,x,x,z,z,z,z) + sympy.integrate(sympy.diff(sphysx, x,x)*sympy.diff(sphysz, z,z) ,x,x,x,x,z,z,z,z)
#test_forward(A, sphys, ssol, xg, 4)

#
# Accuracy tests: qid
#
print("qid:")
A = c3d.qid(nx,ny,nz, 3,3,3)
#sphysx = np.sum([np.random.ranf()*x**i for i in np.arange(0,nx,1)])
#sphysz = np.sum([np.random.ranf()*z**i for i in np.arange(0,nz,1)])
#ssol = sphysx*sphysz
#test_forward(A, sphys, ssol, xg, 3)
