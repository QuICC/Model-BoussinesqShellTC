"""Module provides the functions required for any model"""

from __future__ import division
from __future__ import unicode_literals

verbose_write_mtx = False
if verbose_write_mtx:
    import scipy.io as io

import geomhdiscc.base.utils as utils


class BaseModel:
    """Base class for all the models"""

    linearize = False

    use_galerkin = False

    SOLVER_HAS_BC = 0
    SOLVER_NO_TAU = 1
    STENCIL = 2
    FIELD_TO_RHS = 3

    SLOWEST = 0
    MODE = 1
    SINGLE = 2
    GEOMETRIC_1D_3D = 3

    def time(self, res, eq_params, eigs, bcs, fields):
        """Create the time derivative operator"""

        mat = utils.build_diag_matrix(fields, self.time_block, (res,eq_params,eigs,bcs))
        if verbose_write_mtx:
            fname = "matrix_time_" + str(bcs["bcType"])
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def implicit_linear(self, res, eq_params, eigs, bcs, fields):
        """Create the implicit linear operator"""

        mat = utils.build_block_matrix(fields, self.linear_block, (res,eq_params,eigs,bcs))
        if verbose_write_mtx:
            fname = "matrix_linear_" + str(bcs["bcType"])
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname  + ".mtx", mat)
        return mat

    def explicit_linear(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create the explicit linear operator"""

        mat = -self.linear_block(res, eq_params, eigs, bcs, field_row, field_col)
        if verbose_write_mtx:
            fname = "matrix_explicit_" + str(bcs["bcType"])
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat
