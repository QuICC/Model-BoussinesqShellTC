"""Module provides generic functions for the sparse chebyshev representation."""

from __future__ import division
from __future__ import unicode_literals

import scipy.sparse as spsp


def build_diagonals(ns, nzrow, ds, offsets, cross_parity = None, has_wrap = True):
    """Build diagonals from function list and offsets"""

    # Build wrapped around values
    if has_wrap:
        wrap = [[0]*abs(offsets[0]) for i in range(len(ds))]
        for d,f in enumerate(ds):
            lb = max(0, -offsets[d])
            for i in range(0, lb):
                if ns[i] > nzrow:
                    lbw = sum(offsets < 0)
                    step = (ns[1]-ns[0])
                    if cross_parity == None:
                        shift = 0
                    else:
                        shift = (-1)**cross_parity 
                    col = sum(((ns + shift) - (-(ns[i] + shift) - step*offsets[d])) < 0)
                    col = lbw + col - i
                    row = i - max(0, lbw - col)
                    wrap[col][row] = f(ns[i])

    # Build diagonals
    diags = [0]*len(ds)
    for d,f in enumerate(ds):
        lb = max(0, -offsets[d])
        ub = min(len(ns),len(ns)-offsets[d])
        diags[d] = [f(n) if n > nzrow else 0 for n in ns[lb:ub]]
        if has_wrap:
            if len(wrap[d]) > 0:
                diags[d][0:len(wrap[d])] = [x + y for x, y in zip(diags[d][0:len(wrap[d])], wrap[d][:])]

    return diags

def build_block_matrix(fields, func, func_args):

    tmp = []
    for field_row in fields:
        row = []
        for field_col in fields:
            args = func_args + (field_row,field_col)
            row.append(func(*args))
        tmp.append(row)

    return spsp.bmat(tmp)

def build_diag_matrix(fields, func, func_args):

    tmp = []
    for field_row in fields:
        args = func_args + (field_row,)
        tmp.append(func(*args))
   
    return spsp.block_diag(tmp)

def triplets(mat):
    mat = mat.tocoo();

    return list(zip(mat.row,mat.col,mat.data))

def restricted_kron_2d(A, B, restriction = None):
    """Compute a double Kronecker product with possible restrictions"""

    if restriction == None or A.nnz == 0 or B.nnz == 0:
        mat = spsp.kron(A, B)

    else:
        diag = lil_matrix((1,A.shape[0]))
        diag[0,restriction] = 1.0
        S = spsp.diags(diag.todense(), [0], shape = A.shape)

        mat = spsp.kron(S*A, B)

    return mat

def restricted_kron_3d(A, B, C, restriction = None):
    """Compute a triple Kronecker product with possible restrictions"""

    if restriction == None or A.nnz == 0 or B.nnz == 0 or C.nnz == 0:
        mat = spsp.kron(A, spsp.kron(B, C))

    else:
        output_shape = (B.shape[0]*C.shape[0], A.shape[1]*B.shape[1]*C.shape[1])
        itSlow = iter(restriction[0])
        itFast = iter(restriction[1])
        row = itSlow.next()
        lines = itFast.next()
        if row == 0:
            mat = spsp.kron(A[i,:], restricted_kron_2d(B, C, lines))
            row = itSlow.next()
            lines = itFast.next()
        else:
            mat = spsp.coo_matrix(output_shape)

        for i in range(1, AA.shape[0]):
            if i == row:
                mat = spsp.vstack([mat, spsp.kron(A[i,:], restricted_kron_2d(B, C, lines))])
                row = itSlow.next()
                lines = itFast.next()
            else:
                mat = spsp.vstack([mat, spsp.coo_matrix(output_shape)])

    return mat
