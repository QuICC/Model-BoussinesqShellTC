"""Module provides the functions to generate the Boussinesq convection in a rotating cylindrical annulus model"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.annulus as annulus
from geomhdiscc.geometry.cylindrical.annulus_boundary import no_bc
import geomhdiscc.base.base_model as base_model


class BoussinesqRotConvAnnulus(base_model.BaseModel):
    """Class to setup the Boussinesq convection in a rotating cylindrical annulus model"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["taylor", "prandtl", "rayleigh", "ro", "rratio"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocityx", "velocityy", "velocityz", "pressure", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocitx",""), ("velocityy",""), ("velocityz",""), ("pressure",""), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocityx",""), ("velocityy",""), ("velocityz",""), ("pressure",""), ("temperature","")]

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

        # Impose no boundary conditions
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
        m = eigs[0]

        if field_row == ("velocityx",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, no_bc)

        elif field_row == ("velocityy",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, no_bc)

        elif field_row == ("velocityy",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, no_bc)

        elif field_row == ("temperature",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, no_bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        a, b = self.linear_r2x(eq_params['ro'], eq_params['rratio'])
        m = eigs[0]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = annulus.i2j2x2lapl(res[0], res[1], m, a, b, bc) + annulus.i2j2(res[0], res[1], a, b, bc, -1.0)

            elif field_col == ("velocityy",""):
                mat = annulus.i2j2(res[0], res[1], a, b, bc, -2.0*1j*m)

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x2d1(res[0], res[1], bc)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[1], bc)

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = annulus.i2j2(res[0], res[1], a, b, bc, -2.0*1j*m)

            elif field_col == ("velocityy",""):
                mat = annulus.i2j2x2lapl(res[0], res[1], m, a, b, bc) + annulus.i2j2(res[0], res[1], a, b, bc, -1.0)

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x1(res[0], res[1], bc, 1j*m)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[1], bc)

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("velocityy",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("velocityz",""):
                mat = annulus.i2j2x2lapl(res[0], res[1], m, a, b, bc)

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x2e1(res[0], res[1], bc)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[1], bc)

        elif field_row == ("pressure",""):
            if field_col == ("velocityx",""):
                mat = annulus.i1j1(res[0], res[1], bc) + annulus.i1j1x1d1(res[0], res[1], bc)

            elif field_col == ("velocityy",""):
                mat = annulus.ij1(res[0], res[1],  bc, 1j*m)

            elif field_col == ("velocityz",""):
                mat = annulus.i1j1x1e1(res[0], res[1], bc)

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("pressure",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[1], bc)

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("velocityy",""):
                mat = annulus.zblk(res[0], res[1],  bc)

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("pressure",""):
                mat = annulus.zblk(res[0], res[1], bc)

            elif field_col == ("temperature",""):
                mat = annulus.i2j2x2lapl(res[0], res[1], m, a, b, bc)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Ta = eq_params['taylor']
        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        a, b = self.linear_r2x(eq_params['ro'], eq_params['rratio'])
        m = eigs[0]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, bc)

        elif field_row == ("velocityy",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, bc)

        elif field_row == ("velocityz",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, bc)

        elif field_row == ("pressure",""):
            mat = annulus.zblk(res[0], res[1], bc)

        elif field_row == ("temperature",""):
            mat = annulus.i2j2x2(res[0], res[1], a, b, bc, Pr)

        return mat

    def linear_r2x(self, ro, rratio):
        """Calculat a and b for linear map r = a*x + b"""

        b = (ro*rratio + ro)/2.0;
        a = ro - b;

        return (a, b)
