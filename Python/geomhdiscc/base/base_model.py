"""Module provides the functions required for any model"""

from __future__ import division
from __future__ import unicode_literals

import geomhdiscc.base.utils as utils
#import scipy.io as io

class BaseModel:
    """Base class for all the models"""

    linearize = False

    def time(self, res, eq_params, eigs, bcs, fields):
        """Create the time derivative operator"""

        mat = utils.build_diag_matrix(fields, self.time_block, (res,eq_params,eigs,bcs))
        #io.mmwrite("matrix_time_" + str(bcs["bcType"]) + "_"+ str(eigs[0])  + "_"+ str(eigs[1]) + ".mtx", mat)
        return mat


    def implicit_linear(self, res, eq_params, eigs, bcs, fields):
        """Create the implicit linear operator"""

        mat = utils.build_block_matrix(fields, self.linear_block, (res,eq_params,eigs,bcs))
        #io.mmwrite("matrix_linear_" + str(bcs["bcType"]) + "_"+ str(eigs[0]) + "_"+ str(eigs[1]) + ".mtx", mat)
        return mat


    def explicit_linear(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create the explicit linear operator"""

        mat = -self.linear_block(res, eq_params, eigs, bcs, field_row, field_col)
        #io.mmwrite("matrix_explicit_" + str(bcs["bcType"]) + "_"+ str(eigs[0]) + "_"+ str(eigs[1]) + ".mtx", mat)
        return mat
