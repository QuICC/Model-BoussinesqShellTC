"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a 3D box (streamfunction formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_3d as c3d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_3d import no_bc


class BoussinesqRB3DBoxST(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a 3D box (streamfunction formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "scale1d", "scale2d", "scale3d"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["streamfunction", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("streamfunction",""), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("streamfunction",""), ("temperature","")]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]*res[1]*res[2]
        if self.use_galerkin:
            if field_row == ("streamfunction","") or field_row == ("temperature",""):
                shift_x = 2
                shift_y = 2
                shift_z = 2
            else:
                shift_x = 0
                shift_y = 0
                shift_z = 0

            gal_n = (res[0] - shift_x)*(res[1] - shift_y)*(res[2] - shift_z)

        else:
            gal_n = tau_n
            shift_x = 0
            shift_y = 0
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_x,shift_y,shift_z), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = False

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode:
        index_mode = self.SINGLE

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
            # No-slip/No-slip/No-slip, Fixed temperature/Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("streamfunction",""):
                        bc = {'x':{0:-40, 'r':0}, 'y':{0:-40, 'r':0}, 'z':{0:-40, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("streamfunction","") and field_col == ("streamfunction",""):
                        bc = {'x':{0:40}, 'y':{0:40}, 'z':{0:40}, 'priority':'xz'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'zy'}

            # Stress-free/Stress-free/no-slip, Fixed flux/Fixed flux/Fixed temperature
            elif bcId == 4:
                if self.use_galerkin:
                    if field_col == ("streamfunction",""):
                        bc = {'x':{0:-41, 'r':0}, 'y':{0:-41, 'r':0}, 'z':{0:-40, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("streamfunction","") and field_col == ("streamfunction",""):
                        bc = {'x':{0:41}, 'y':{0:41}, 'z':{0:40}, 'priority':'zsx'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'y':{0:21}, 'z':{0:20}, 'priority':'zsx'}

            # Stress-free/Stress-free/no-slip, Fixed flux/Fixed flux/Fixed temperature
            elif bcId == 6:
                if self.use_galerkin:
                    if field_col == ("streamfunction",""):
                        bc = {'x':{0:-41, 'r':0}, 'y':{0:-41, 'r':0}, 'z':{0:-40, 'r':0}}

                else:
                    if field_row == ("streamfunction","") and field_col == ("streamfunction",""):
                        bc = {'x':{0:41}, 'y':{0:41}, 'z':{0:41}, 'priority':'zsx'}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("streamfunction",""):
                    bc['x']['r'] = 4
                    bc['y']['r'] = 4
                    bc['z']['r'] = 4
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("streamfunction",""):
                        bc = {'x':{0:-40, 'r':0}, 'y':{0:-40, 'r':0}, 'z':{0:-40, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 4:
                    if field_col == ("streamfunction",""):
                        bc = {'x':{0:-41, 'r':0}, 'y':{0:-41, 'r':0}, 'z':{0:-40, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("streamfunction",""):
                    bc['x']['r'] = 4
                    bc['y']['r'] = 4
                    bc['z']['r'] = 4
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return c3d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create the quasi-inverse operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("streamfunction",""):
            mat = c3d.i4j4k4(res[0], res[1], res[2], bc)

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Ra = eq_params['rayleigh']

        xscale = eq_params['scale1d']
        yscale = eq_params['scale2d']
        zscale = eq_params['scale3d']

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("streamfunction",""):
            if field_col == ("streamfunction",""):
                mat = c3d.i4j4k4lapl2(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale)

            elif field_col == ("temperature",""):
                mat = c3d.i4j4k4(res[0], res[1], res[2], bc, -Ra/16.0)

        elif field_row == ("temperature",""):
            if field_col == ("streamfunction",""):
                mat = c3d.i2j2k2laplh(res[0], res[1], res[2], bc, -1.0, xscale = xscale, yscale = yscale)

            elif field_col == ("temperature",""):
                mat = c3d.i2j2k2lapl(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']

        xscale = eq_params['scale1d']
        yscale = eq_params['scale2d']
        zscale = eq_params['scale3d']

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("streamfunction",""):
            mat = c3d.i4j4k4lapl(res[0], res[1], res[2], bc, 1.0/Pr, xscale = xscale, yscale = yscale, zscale = zscale)

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc)

        return mat
