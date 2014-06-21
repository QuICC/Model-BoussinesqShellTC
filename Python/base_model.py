"""Module provides the functions required for any model"""

from __future__ import division

import utils

class BaseModel:
   """Base class for all the models"""

   def time(self, res, eq_params, eigs, bcs, fields):
      """Create the time derivative operator"""

      mat = utils.build_diag_matrix(fields, self.time_block, (res,eq_params,eigs,bcs))
      return mat


   def implicit_linear(self, res, eq_params, eigs, bcs, fields):
      """Create the implicit linear operator"""

      mat = utils.build_block_matrix(fields, self.linear_block, (res,eq_params,eigs,bcs))
      return mat


   def explicit_linear(self, res, eq_params, eigs, bcs, field_row, field_col):
      """Create the explicit linear operator"""

      mat = linear_block(res, eq_params, eigs, self.field_row, field_col)
      return mat
