"""Module provides the functions to generate the test model for the CFT scheme"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.cylinder as cylinder
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cylindrical.cylinder_boundary import no_bc


class TestCFTScheme(base_model.BaseModel):
    """Class to setup the test model for the CFT scheme"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocityx", "velocityy", "velocityz", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields = [("velocityx",""), ("velocityy",""), ("velocityz",""), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        # Solve as coupled equations
        fields = [("velocityx",""), ("velocityy",""), ("velocityz",""), ("temperature","")]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        return []

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]*res[2]
        if self.use_galerkin:
            if field_row == ("velocityx",""):
                shift_r = 2
                shift_z = 2
            elif field_row == ("velocityy",""):
                shift_r = 2
                shift_z = 2
            elif field_row == ("velocityz",""):
                shift_r = 2
                shift_z = 2
            elif field_row == ("temperature",""):
                shift_r = 2
                shift_z = 2
            else:
                shift_r = 0
                shift_z = 0

            gal_n = (res[0] - shift_r)*(res[2] - shift_z)

        else:
            gal_n = tau_n
            shift_r = 0
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_r,0,shift_z), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is real
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: SLOWEST, MODE, GEOMETRIC_1D_3D
        index_mode = self.SLOWEST

        # Compute block info
        block_info = self.block_size(res, field_row)

        # Compute system size
        sys_n = 0
        for f in im_fields:
            sys_n += self.block_size(res, f)[1]
        
        if sys_n == 0:
            sys_n = block_info[1]
        block_info = block_info + (sys_n,)

        return (is_complex, im_fields, ex_fields, index_mode, block_info)

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
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'r':{0:10}, 'z':{0:20}}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'r':{0:10}, 'z':{0:20}}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'r':{0:10}, 'z':{0:20}}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'r':{0:10}, 'z':{0:20}}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("velocityy",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityy",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
        
        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("velocityz",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2

        else:
            bc = no_bc()

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return cylinder.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        m = eigs[0]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            mat = cylinder.i2j2x2(res[0],res[2], (m+1)%2, bc)

        elif field_row == ("velocityy",""):
            mat = cylinder.i2j2x2(res[0],res[2], (m+1)%2, bc)

        elif field_row == ("velocityz",""):
            mat = cylinder.i2j2x2(res[0],res[2], m%2, bc)

        elif field_row == ("temperature",""):
            mat = cylinder.i2j2x2(res[0],res[2], m%2, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block of linear operator"""

        m = eigs[0]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = cylinder.i2j2x2lapl(res[0], res[2], m, (m+1)%2, bc)
                bc['r'][0] = min(bc['r'][0], 0)
                bc['z'][0] = min(bc['z'][0], 0)
                mat = mat + cylinder.i2j2(res[0], res[2], (m+1)%2, bc, -1.0)

            elif field_col == ("velocityy",""):
                mat = cylinder.i2j2(res[0], res[2], (m+1)%2, bc, -2.0*1j*m)

            elif field_col == ("velocityz",""):
                mat = cylinder.zblk(res[0], res[2], m%2, 1, 2, bc)

            elif field_col == ("temperature",""):
                mat = cylinder.zblk(res[0], res[2], m%2, 1, 2, bc)

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = cylinder.i2j2(res[0], res[2], (m+1)%2, bc, 2.0*1j*m)

            elif field_col == ("velocityy",""):
                mat = cylinder.i2j2x2lapl(res[0], res[2], m, (m+1)%2, bc)
                bc['r'][0] = min(bc['r'][0], 0)
                bc['z'][0] = min(bc['z'][0], 0)
                mat = mat + cylinder.i2j2(res[0], res[2], (m+1)%2, bc, -1.0)

            elif field_col == ("velocityz",""):
                mat = cylinder.zblk(res[0], res[2], m%2, 1, 2, bc)

            elif field_col == ("temperature",""):
                mat = cylinder.zblk(res[0], res[2], m%2, 1, 2, bc)

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = cylinder.zblk(res[0], res[2], (m+1)%2, 1, 2, bc)

            elif field_col == ("velocityy",""):
                mat = cylinder.zblk(res[0], res[2], (m+1)%2, 1, 2, bc)

            elif field_col == ("velocityz",""):
                mat = cylinder.i2j2x2lapl(res[0], res[2], m, m%2, bc)

            elif field_col == ("temperature",""):
                mat = cylinder.zblk(res[0], res[2], m%2, 1, 2, bc)

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = cylinder.zblk(res[0], res[2], (m+1)%2, 1, 2, bc)

            elif field_col == ("velocityy",""):
                mat = cylinder.zblk(res[0], res[2], (m+1)%2, 1, 2, bc)

            elif field_col == ("velocityz",""):
                mat = cylinder.zblk(res[0], res[2],  m%2, 1, 2, bc)

            elif field_col == ("temperature",""):
                mat = cylinder.i2j2x2lapl(res[0], res[2], m, m%2, bc)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        m = eigs[0]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = cylinder.i2j2x2(res[0],res[2], (m+1)%2, bc)

        elif field_row == ("velocityy",""):
            mat = cylinder.i2j2x2(res[0],res[2], (m+1)%2, bc)

        elif field_row == ("velocityz",""):
            mat = cylinder.i2j2x2(res[0],res[2], m%2, bc)

        elif field_row == ("temperature",""):
            mat = cylinder.i2j2(res[0],res[2], m%2, bc)

        return mat