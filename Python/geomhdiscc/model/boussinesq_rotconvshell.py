"""Module provides the functions to generate the Boussinesq convection in a rotating spherical shell model"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.shell as sphere
import geomhdiscc.base.base_model as base_model


class BoussinesqRotConvShell(base_model.BaseModel):
    """Class to setup the Boussinesq convection in a rotating spherical shell model"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["ekman", "prandtl", "rayleigh", "ro", "rratio"]

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

        fields =  [("velocity","tor"), ("velocity","pol"), ("temperature","")]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Equation doesn't have geometric coupling
        has_geometric_coupling = True
        # Index mode: SLOWEST = 0, MODE = 1
        index_mode = 1

        # Rows per equation block and number of rhs
        block_info = (res[0]*res[1], 1)

        return (is_complex, im_fields, ex_fields, has_geometric_coupling, index_mode, block_info)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        use_tau_boundary = True

        # Impose no boundary conditions
        no_bc = [0]
        if bcs["bcType"] == 2:
            bc = no_bc
        else:
            # Impose no boundary conditions
            if bcs["bcType"] == 1 and use_tau_boundary:
                bc = no_bc
            else: #bcType == 0 or Galerkin boundary
                l = eigs[0]
                m = eigs[1]

                bc = None
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    bc_field = {}
                    bc_field[("velocity","tor")] = [20]
                    bc_field[("velocity","pol")] = [40]
                    bc_field[("temperature","")] = [20]

                    if field_col == field_row:
                        bc = bc_field[field_col]
                elif bcId == 1:
                    bc_field = {}
                    bc_field[("velocity","tor")] = [21]
                    bc_field[("velocity","pol")] = [41]

                    if field_col == field_row:
                        bc = bc_field[field_col]

                if bc is None:
                    if use_tau_boundary:
                        bc = no_bc
                    else:
                        bc = bc_field[field_col]
                        bc[0] = -bc[0]

        return bc

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        a, b = self.linear_r2x(eq_params['ro'], eq_params['rratio'])
        l = eigs[0]
        m = eigs[1]

        if field_row == ("velocity","tor"):
            mat = shell.i2x2(res[0], res[1], m, a, b, [0])

        elif field_row == ("velocity","pol"):
            mat = shell.i4x4(res[0], res[1], m, a, b, [0])

        elif field_row == ("temperature",""):
            mat = shell.i2x2(res[0], res[1], m, a, b, [0])

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        a, b = self.linear_r2x(eq_params['ro'], eq_params['rratio'])
        l = eigs[0]
        m = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = shell.i2x2lapl(res[0], res[1], m, a, b, bc, 1.0, 'laplh') + shell.i2x2(res[0], res[1], m, a, b, [min(bc[0],0)], 1j*m*Ta**0.5)

            elif field_col == ("velocity","pol"):
                mat = shell.i2x2coriolis(res[0], res[1], m, a, b, bc, -Ta**0.5)

            elif field_col == ("temperature",""):
                mat = shell.zblk(res[0], res[1], m, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = shell.i4x4coriolis(res[0], res[1], m, a, b, bc, Ta**0.5)

            elif field_col == ("velocity","pol"):
                mat = shell.i4x4lapl2(res[0], res[1], m, a, b, bc, 1.0, 'laplh') + shell.i4x4lapl(res[0], res[1], m, a, b, [min(bc[0],0)], 1j*m*Ta**0.5)

            elif field_col == ("temperature",""):
                mat = shell.i4x4(res[0], res[1], m, a, b, bc, -Ra, 'laplh')

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = shell.zblk(res[0], res[1], m, bc)

            elif field_col == ("velocity","pol"):
                if self.linearize:
                    mat = shell.i2x2(res[0], res[1], m, a, b, bc, 1.0, 'laplh')

                else:
                    mat = shell.zblk(res[0], res[1], m, bc)

            elif field_col == ("temperature",""):
                mat = shell.i2x2lapl(res[0], res[1], m, a, b, bc)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        a, b = self.linear_r2x(eq_params['ro'], eq_params['rratio'])
        l = eigs[0]
        m = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = shell.i2x2(res[0], res[1], m, a, b, bc, 1.0, 'laplh')

        elif field_row == ("velocity","pol"):
            mat = shell.i4x4lapl(res[0], res[1], m, a, b, bc, 1.0, 'laplh')

        elif field_row == ("temperature",""):
            mat = shell.i2x2(res[0], res[1], m, a, b, bc, Pr)

        return mat

    def linear_r2x(self, ro, rratio):
        """Calculat a and b for linear map r = a*x + b"""

        b = (ro*rratio + ro)/2.0;
        a = ro - b;

        return (a, b)
