"""Module provides the functions to generate the test model for the TFT scheme"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_2d as geo
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_2d import no_bc


class TestTFTScheme(base_model.BaseModel):
    """Class to setup the test model for the TFT scheme"""

    solve_coupled = True

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "scale1d", "scale3d"]

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields = [("velocity","x"), ("velocity","y"), ("velocity","z"), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        # Solve as coupled equations
        if self.solve_coupled:
            fields = [("velocity","x"), ("velocity","y"), ("velocity","z"), ("temperature","")]

        # Solve as splitted equations
        else:
            fields = [field_row]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row in [("velocity","x"), ("velocity","y"), ("velocity","z"), ("temperature","")]:
                fields = [field_row]
            else:
                fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]*res[2]
        if self.use_galerkin:
            if field_row == ("velocityx",""):
                shift_x = 2
                shift_z = 2
            elif field_row == ("velocityy",""):
                shift_x = 2
                shift_z = 2
            elif field_row == ("velocityz",""):
                shift_x = 2
                shift_z = 2
            elif field_row == ("temperature",""):
                shift_x = 2
                shift_z = 2
            else:
                shift_x = 0
                shift_z = 0

            gal_n = (res[0] - shift_x)*(res[2] - shift_z)

        else:
            gal_n = tau_n
            shift_x = 0
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_x,0,shift_z), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is real
        is_complex = True

        # Index mode: 
        index_mode = self.SLOWEST_SINGLE_RHS

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:20}}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:20}, 'z':{0:20}}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:20}, 'z':{0:20}}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'z':{0:20}}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("velocityy",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'rt':0}, 'z':{0:-20, 'rt':0}}
        
        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['rt'] = 2
                    bc['z']['rt'] = 2

        else:
            bc = no_bc()

        return bc

    def nonlinera_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create the explicit nonlinear operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","x") and field_col == field_row:
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("velocity","y") and field_col == field_row:
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("velocity","z") and field_col == field_row:
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("temperature","") and field_col == field_row:
            mat = geo.i2j2(res[0],res[2], bc)

        else:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block of linear operator"""

        k = eigs[0]/2.0
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = geo.i2j2lapl(res[0], res[2], k, bc)

            elif field_col == ("velocityy",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = geo.i2j2lapl(res[0], res[2], k, bc)

            elif field_col == ("velocityz",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = geo.i2j2lapl(res[0], res[2], k, bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = geo.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = geo.i2j2lapl(res[0], res[2], k, bc)

        else:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("velocityy",""):
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("velocityz",""):
            mat = geo.i2j2(res[0],res[2], bc)

        elif field_row == ("temperature",""):
            mat = geo.i2j2(res[0],res[2], bc)

        else:
            raise RuntimeError("Equations are not setup properly!")

        return mat
