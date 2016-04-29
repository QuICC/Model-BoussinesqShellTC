"""Module provides functions to generate the boundary conditions in a cylinder with Worland expansion in radius"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.geometry.cartesian.cartesian_boundary_1d as c1dbc
import geomhdiscc.geometry.cylindrical.cylinder_radius_boundary_worland as radbc


def no_bc():
    """Get a no boundary condition flag"""

    return {'r':radbc.no_bc(), 'z':c1dbc.no_bc()}

def brid(n, m, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(n-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = radbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return radbc.constrain(mat, m, tbc)

def bzid(n, q, d, bc, location = 't'):
    """Create a boundary indentity"""

    if bc[0] < 0:
        mat = spsp.eye(n-q, n-(-bc[0])//10)
    else:
        offsets = [d]
        diags = [[1]*(n-abs(d))]

        mat = spsp.diags(diags, offsets).tolil()
        if location == 't':
            mat[0:q,:] = 0
        elif location == 'b':
            if q > 0:
                mat[-q:,:] = 0

    tbc = c1dbc.no_bc()
    for key, val in bc.items():
        if key != 0:
            tbc[key] = val

    return c1dbc.constrain(mat,tbc)

def constrain(mat, nr, nz, m, qr, qz, bc, location = 't'):
    """Contrain the matrix with the Tau boundary condition"""

    priority = bc.get('priority', 'r')
    sr = 0
    dr = 0
    sz = 0
    dz = 0
    if priority == 'r':
        sr = qr
    elif priority == 'z':
        sz = qz
    elif priority == 'n':
        sr = qr
        sz = qz
    elif priority == 'sr':
        sr = qr
        dr = -qr
    elif priority == 'sz':
        sz = qz
        dz = -qz
    else:
        raise RuntimeError("Unknown boundary condition priority!")

    bc_mat = mat
    if bc['r'][0] > 0:
        bcMat = spsp.lil_matrix((nr,nr))
        bcMat = radbc.constrain(bcMat, m, bc['r'], location = location)
        bc_mat = bc_mat + spsp.kron(bzid(nz, sz, 0, bc['z'], location = location), bcMat, format = 'coo')

    if bc['z'][0] > 0:
        bcMat = spsp.lil_matrix((nz,nz))
        bcMat = c1dbc.constrain(bcMat, bc['z'], location = location)
        bc_mat = bc_mat + spsp.kron(bcMat, brid(nr, m, sr, 0, bc['r'], location = location), format = 'coo')

    return bc_mat
