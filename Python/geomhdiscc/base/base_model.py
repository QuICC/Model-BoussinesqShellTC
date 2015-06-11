"""Module provides the functions required for any model"""

from __future__ import division
from __future__ import unicode_literals

verbose_write_mtx = False
if verbose_write_mtx:
    import scipy.io as io
    import os

import geomhdiscc.base.utils as utils


class BaseModel:
    """Base class for all the models"""

    linearize = False

    use_galerkin = False

    SOLVER_HAS_BC = 0
    SOLVER_NO_TAU = 1
    STENCIL = 2
    FIELD_TO_RHS = 3

    SLOWEST_SINGLE_RHS = 0
    SLOWEST_MULTI_RHS = 1
    MODE = 2
    SINGLE = 3

    EXPLICIT_LINEAR = 0
    EXPLICIT_NONLINEAR = 1
    EXPLICIT_NEXTSTEP = 2

    def time(self, res, eq_params, eigs, bcs, fields, restriction = None):
        """Create the time derivative operator"""

        mat = utils.build_diag_matrix(fields, self.time_block, (res,eq_params,eigs,bcs), restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_time_" + str(bcs["bcType"]) + "_" + str(os.getpid())
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def implicit_linear(self, res, eq_params, eigs, bcs, fields, restriction = None):
        """Create the implicit linear operator"""

        mat = utils.build_block_matrix(fields, self.implicit_block, (res,eq_params,eigs,bcs), restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_linear_" + str(bcs["bcType"]) + "_" + str(os.getpid())
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname  + ".mtx", mat)
        return mat

    def explicit_linear(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.explicit_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_linear_" + ''.join(field_row) + '_' + ''.join(field_col) + '_' + str(bcs["bcType"]) + "_" + str(os.getpid())
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def explicit_nonlinear(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.nonlinear_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_nonlinear_" + ''.join(field_row) + '_' + ''.join(field_col) + '_' + str(bcs["bcType"]) + "_" + str(os.getpid())
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def explicit_nextstep(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit linear operator"""

        mat = self.nextstep_block(res, eq_params, eigs, bcs, field_row, field_col, restriction = restriction)
        if verbose_write_mtx:
            fname = "matrix_explicit_nextstep_" + ''.join(field_row) + '_' + ''.join(field_col) + '_' + str(bcs["bcType"]) + "_" + str(os.getpid())
            for e in eigs:
                fname = fname + "_" + str(e)
            io.mmwrite(fname + ".mtx", mat)
        return mat

    def compile_equation_info(self, res, field_row, is_complex, index_mode):
        """Collect all equation info together"""

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        lin_fields = self.explicit_fields(self.EXPLICIT_LINEAR, field_row)
        # Additional explicit nonlinear fields
        nl_fields = self.explicit_fields(self.EXPLICIT_NONLINEAR, field_row)
        # Additional explicit update for next step linear fields
        next_fields = self.explicit_fields(self.EXPLICIT_NEXTSTEP, field_row)

        # Compute block info
        block_info = self.block_size(res, field_row)

        # Compute system size
        sys_n = 0
        for f in im_fields:
            sys_n += self.block_size(res, f)[1]
        
        if sys_n == 0:
            sys_n = block_info[1]
        block_info = block_info + (sys_n,)

        return (is_complex, im_fields, lin_fields, nl_fields, next_fields, index_mode, block_info)

    def stability_sizes(self, res, eigs):
        """Get the block sizes in the stability calculation matrix"""

        # Block sizes
        blocks = []
        for f in self.stability_fields():
            blocks.append(self.block_size(res, f)[1])

        # Invariant size (local dimension in spectral space, no restriction)
        invariant = (res[0],)*len(self.stability_fields())

        # Index shift
        shift = 0

        return (blocks, invariant, shift)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        raise NotImplementedError("Model cannot be used for linear stability calculations!")

    def block_size(self, res, field_row):
        """Create block size information"""

        raise NotImplementedError("Model cannot be used for linear stability calculations!")

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        raise NotImplementedError("Model cannot be used for linear stability calculations!")

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        raise NotImplementedError("Model should implement this method!")

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        raise NotImplementedError("Model should implement this method!")

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for implicit operator"""

        raise NotImplementedError("Model should implement this method!")

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        raise NotImplementedError("Model should implement this method!")

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nonlinear term"""

        raise NotImplementedError("Model should implement this method!")

    def nextstep_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nextstep update"""

        raise NotImplementedError("Model should implement this method!")

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""
        
        raise NotImplementedError("Stencil needs to be implemented in model!")
