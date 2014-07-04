"""Module provides the functions to generate the test model for the TTT scheme"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
from geomhdiscc.base.utils import triplets
import geomhdiscc.geometry.cartesian.cartesian_3d as c3d
import geohmdiscc.base.base_model as base_model


class TestTTT(base_model.BaseModel):
    """Class to setup the test model for the TTT scheme"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "gamma", "chi"]


    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]


    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["streamfunction", "velocityz", "temperature"]


    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        # Solve as coupled equations
        if True:
            fields = [("streamfunction",""), ("velocityz",""), ("temperature","")]

        # Solve as splitted equations
        else:
            fields = [field_row]

        return fields


    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        return []


    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is real
        is_complex = False

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Equation doesn't have geometric coupling
        has_geometric_coupling = False
        # Index mode: SLOWEST = 0, MODE = 1
        index_mode = 0

        # Rows per equation block and number of rhs
        block_info = (res[0]*res[1]*res[2], 1)

        return (is_complex, im_fields, ex_fields, has_geometric_coupling, index_mode,block_info)


    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        use_tau_boundary = True
        # Impose no boundary conditions
        no_bc = {'x':[0],'y':[0],'z':[0]}
        if bcs["bcType"] == 2:
            bc = no_bc
        else:
            # Impose no boundary conditions
            if bcs["bcType"] == 1 and use_tau_boundary:
                bc = no_bc
            else: #bcType == 0 or Galerkin boundary
                bc = None
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    bc_field = {}
                    bc_field[("streamfunction","")] = {'x':[40],'y':[40],'z':[40]}
                    bc_field[("velocityz","")] = {'x':[20],'y':[20],'z':[20]}
                    bc_field[("temperature","")] = {'x':[20],'y':[20],'z':[20]}
                    if field_col == field_row:
                        bc = bc_field[field_col]

                if bc is None:
                    if use_tau_boundary:
                        bc = no_bc
                    else:
                        bc = {}
                        for k,v in bc_field[field_col]:
                            bc[k] = v
                            bc[k][0] = -v[0]

        return bc


    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        if field_row == ("streamfunction",""):
            mat = c3d.i2j2k2(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

        elif field_row == ("velocityz",""):
            mat = c3d.i4j4k4(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0],res[1],res[2], {'x':[0], 'y':[0], 'z':[0]})

        return mat


    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block of linear operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("streamfunction",""):
            if field_row == ("streamfunction",""):
                mat = c3d.i2j2k2lapl(res[0],res[1],res[2], bc)

            elif field_row == ("velocityz",""):
                mat = c3d.zblk(res[0],res[1],res[2],4,4,4, bc)

            elif field_row == ("temperature",""):
                mat = c3d.zblk(res[0],res[1],res[2],4,4,4, bc)

        elif field_row == ("velocityz",""):
            if field_row == ("streamfunction",""):
                mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

            elif field_row == ("velocityz",""):
                mat = c3d.i4j4k4lapl2(res[0],res[1],res[2], bc)

            elif field_row == ("temperature",""):
                mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

        elif field_row == ("temperature",""):
            if field_row == ("streamfunction",""):
                mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

            elif field_row == ("velocityz",""):
                mat = c3d.zblk(res[0],res[1],res[2],2,2,2, bc)

            elif field_row == ("temperature",""):
                mat = c3d.i2j2k2lapl(res[0],res[1],res[2], bc)

        return mat


    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("streamfunction",""):
            mat = c3d.i2j2k2(res[0],res[1],res[2], bc)

        elif field_row == ("velocityz",""):
            mat = c3d.i4j4k4lapl(res[0],res[1],res[2], bc)

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0],res[1],res[2], bc)

        return mat
