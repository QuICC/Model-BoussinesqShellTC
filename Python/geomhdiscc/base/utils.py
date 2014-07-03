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
