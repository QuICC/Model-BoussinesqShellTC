"""Module provides functions to generate the boundary conditions in a cartesian 2D geometry"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import numpy.polynomial.chebyshev as cheby
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import geomhdiscc.transform.cartesian as phys


def no_bc():
    """Get a no boundary condition flag"""

    return {'x':c1dbc.no_bc(), 'z':c1dbc.no_bc()}

def bid(nx, q, d, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [-d]
        diags = [[1]*(nx-d)]

        mat = spsp.diags(diags, offsets).tolil()
        mat[0:q,:] = 0

    tbc = c1dbc.no_bc()
    tbc['cr'] = bc.get('cr', 0)
    tbc['rb'] = bc.get('rb', 0)
    tbc['cl'] = bc.get('cl', 0)
    tbc['rt'] = bc.get('rt', 0)
    return c1dbc.constrain(mat, tbc)

def sid(nx, q, d, bc):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [0]
        diags = [[1]*(nx-q) + [0]*q]

        mat = spsp.diags(diags, offsets)

    tbc = c1dbc.no_bc()
    tbc['cr'] = bc.get('cr', 0)
    tbc['rb'] = bc.get('rb', 0)
    tbc['cl'] = bc.get('cl', 0)
    tbc['rt'] = bc.get('rt', 0)
    return c1dbc.constrain(mat, tbc)

def bgrid(n, q, d, bc):
    """Create a boundary grid matrix"""

    g = np.cos(np.pi*(np.arange(0,n)+0.5)/n)
    #g = np.cos(np.pi*np.arange(0,n)/(n - 1.0))
    print(g)
    if q > 0:
        g = g[q//2:-q//2]
    print(g)
    mat = np.zeros((n,n))
    for i,x in enumerate(g):
        for j in range(0,n):
            c = np.zeros((n,))
            c[j] = 1
            norm = (1.0 + (j == 0))/2.0
            mat[i+q,j] = norm*cheby.chebval(x,c)

    #mat[0:q,:] = 0

    return  mat

def constrain(mat, nx, nz, qx, qz, bc, location = 't'):
    """Contrain the matrix with the Tau boundary condition"""

    sx = qx
    sz = 0
            
    bc_mat = mat
    if bc['x'][0] > 0:
        bcMat = spsp.lil_matrix((nx,nx))
        bcMat = c1dbc.constrain(bcMat, bc['x'], location = location)
        if location == 't':
            bc_mat = bc_mat + spsp.kron(bid(nz,sz,0,bc['z']), bcMat)
        elif location == 'b':
            bc_mat = bc_mat + spsp.kron(sid(nz,sz,0,bc['z']), bcMat)

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], location = location)
        if bc['x'][0] >= 0:
            if location == 't':
                bc_mat = bc_mat + spsp.kron(bcMat, bid(nx,sx,0,bc['x']))
            elif location == 'b':
                bc_mat = bc_mat + spsp.kron(bcMat, sid(nx,sx,0,bc['x']))
        else:
            if location == 't':
                tmpB = c1dbc.constrain(bid(nx,0,0,c1dbc.no_bc()),bc['x'], location = location)
            elif location == 'b':
                tmpB = c1dbc.constrain(sid(nx,0,0,c1dbc.no_bc()),bc['x'], location = location)
            bc_mat = bc_mat + spsp.kron(bcMat, tmpB)

    return bc_mat
