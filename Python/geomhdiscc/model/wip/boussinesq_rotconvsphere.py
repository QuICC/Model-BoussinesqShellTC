"""Module provides the functions to generate the Boussinesq convection in a rotating sphere model"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.sphere as sphere
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.spherical.sphere_boundary import no_bc


class BoussinesqRotConvSphere(base_model.BaseModel):
    """Class to setup the Boussinesq convection in a rotating sphere model"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["taylor", "prandtl", "rayleigh"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","tor"), ("velocity","pol"), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        if field_row == ("velocity","tor") or field_row == ("velocity","pol") or field_row == ("temperature",""):
            fields =  [("velocity","tor"), ("velocity","pol"), ("temperature","")]
        else:
            fields = []

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]*res[1]
        if self.use_galerkin:
            if field_row == ("velocity","tor") or field_row == ("temperature",""):
                shift_r = 1
            elif field_row == ("velocity","pol"):
                shift_r = 2
            else:
                shift_r = 0

            gal_n = (res[0] - shift_x)*res[1]

        else:
            gal_n = tau_n
            shift_x = 0

        block_info = (tau_n, gal_n, (shift_r,0,0), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: SLOWEST = 0, MODE = 1
        index_mode = self.GEOMETRIC_1D_3D

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
                    if field_col == ("velocity","tor"):
                        bc = {0:-10, 'r':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-20, 'r':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-10, 'r':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:10}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:20}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                            bc = {0:10}

            elif bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-11, 'r':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-21, 'r':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:11}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:21}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['r'] = 1
                elif field_row == ("velocity","pol"):
                    bc['r'] = 2
                elif field_row == ("temperature",""):
                    bc['r'] = 1

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {0:-10, 'r':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-20, 'r':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-10, 'r':0}

                elif bcId == 1:
                    if field_col == ("velocity","tor"):
                        bc = {0:-11, 'r':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-21, 'r':0}
        
        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['r'] = 1
                elif field_row == ("velocity","pol"):
                    bc['r'] = 2
                elif field_row == ("temperature",""):
                    bc['r'] = 1

        else:
            bc = no_bc()

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return shell.stencil(res[0], res[1], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        m = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = sphere.i2x2(res[0], res[1], m, bc)

        elif field_row == ("velocity","pol"):
            mat = sphere.i4x4(res[0], res[1], m, bc)

        elif field_row == ("temperature",""):
            mat = sphere.i2x2(res[0], res[1], m, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        m = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = sphere.i2x2lapl(res[0], res[1], m, bc, 1.0, 'laplh')
                bc[0] = min(bc[0], 0)
                mat = mat + sphere.i2x2(res[0], res[1], m, bc, 1j*m*Ta**0.5)

            elif field_col == ("velocity","pol"):
                mat = sphere.i2x2coriolis(res[0], res[1], m, bc, -Ta**0.5)

            elif field_col == ("temperature",""):
                mat = sphere.zblk(res[0], res[1], m, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = sphere.i4x4coriolis(res[0], res[1], m, bc, Ta**0.5)

            elif field_col == ("velocity","pol"):
                mat = sphere.i4x4lapl2(res[0], res[1], m, bc, 1.0, 'laplh')
                bc[0] = min(bc[0], 0)
                mat = mat + sphere.i4x4lapl(res[0], res[1], m, bc, 1j*m*Ta**0.5)

            elif field_col == ("temperature",""):
                mat = sphere.i4x4(res[0], res[1], m, bc, -Ra, 'laplh')

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = sphere.zblk(res[0], res[1], m, bc)

            elif field_col == ("velocity","pol"):
                if self.linearize:
                    mat = sphere.i2x2(res[0], res[1], m, bc, 1.0, 'laplh')

                else:
                    mat = sphere.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = sphere.i2x2lapl(res[0], res[1], m, bc)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        m = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = sphere.i2x2(res[0], res[1], m, bc, 1.0, 'laplh')

        elif field_row == ("velocity","pol"):
            mat = sphere.i4x4lapl(res[0], res[1], m, bc, 1.0, 'laplh')

        elif field_row == ("temperature",""):
            mat = sphere.i2x2(res[0], res[1], m, bc, Pr)

        return mat