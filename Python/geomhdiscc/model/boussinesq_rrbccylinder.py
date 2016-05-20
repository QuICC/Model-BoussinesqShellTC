"""Module provides the functions to generate the Boussinesq rotating Rayleigh-Benard convection in a cylinder (Toroidal/Poloidal formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp
import functools

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.cylinder_worland as geo
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cylindrical.cylinder_boundary_worland import no_bc


class BoussinesqRRBCCylinder(base_model.BaseModel):
    """Class to setup the Boussinesq rotating Rayleigh-Benard convection in a cylinder (Toroidal/Poloidal formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["taylor", "prandtl", "rayleigh", "gamma"]

    def config_fields(self):
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

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row in [("velocity","tor"), ("velocity","pol"), ("temperature","")]:
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
            if field_row == ("velocity","tor"):
                shift_r = 2
                shift_z = 2
            elif field_row == ("velocity","pol"):
                shift_r = 3
                shift_z = 4
            elif field_row == ("temperature",""):
                shift_r = 1
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

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_SINGLE_RHS

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            m = eigs[0]
            Ra = eq_params['rayleigh']
            zscale = 2

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip/No-slip, Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("temperature",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}

                else:
                    if field_row == ("velocity","tor") and field_col == field_row:
                        bc = {'r':{0:11, 'mixed':{0:15, 'pad':1, 'kron_shift':2, 'kron':functools.partial(geo.c1d.i2d1, cscale=zscale)}}, 'z':{0:20}, 'priority':'r'}
                    elif field_row == ("velocity","tor") and field_col == ("velocity","pol"):
                        bc = {'r':{0:0, 'mixed':{0:16, 'c':-1j*m, 'pad':1, 'kron_shift':2, 'kron':geo.c1d.i2}}, 'z':{0:0}, 'priority':'r'}
                    elif field_row == ("velocity","tor") and field_col == ("temperature",""):
                        bc = {'r':{0:0, 'mixed':{0:10, 'pad':1, 'kron_shift':2, 'kron':geo.c1d.i2, 'c':1j*m*Ra}}, 'z':{0:0}, 'priority':'r'}
                    elif field_row == ("velocity","pol") and field_col == field_row:
                        bc = {'r':{0:22, 'mixed':{0:11, 'pad':2, 'kron_shift':1, 'kron':functools.partial(geo.c1d.i1d1, cscale=zscale)}}, 'z':{0:40}, 'priority':'r'}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","tor"):
                        bc = {'r':{0:0, 'mixed':{0:10, 'c':1j*m*Ra, 'pad':2, 'kron_shift':1, 'kron':geo.c1d.i1}}, 'z':{0:0}, 'priority':'r'}
                    elif field_row == ("temperature","") and field_col == field_row:
                        bc = {'r':{0:10}, 'z':{0:20}, 'priority':'sz'}

            # Stress-free/No-slip, Fixed flux/Fixed temperature
            elif bcId == 2:
                if self.use_galerkin:
                    if field_col == ("temperature",""):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}

                else:
                    if field_row == ("temperature","") and field_col == field_row:
                        bc = {'r':{0:11}, 'z':{0:20}, 'priority':'sz'}
                    elif field_row == ("velocity","tor") and field_col == ("temperature",""):
                        bc = {'r':{0:0, 'mixed':{0:10, 'pad':1, 'c':1j*m*Ra, 'kron_shift':2, 'kron':geo.c1d.i2}}, 'z':{0:0}, 'priority':'r'}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['r']['rt'] = 2
                    bc['z']['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['r']['rt'] = 3
                    bc['z']['rt'] = 4
                elif field_row == ("temperature",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("velocity","pol"):
                        bc = {'r':{0:-10, 'rt':0}, 'z':{0:-20, 'rt':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-11, 'rt':0}, 'z':{0:-20, 'rt':0}}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['r']['rt'] = 1
                    bc['z']['rt'] = 2

        return bc

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        m = eigs[0]

        idx_u, idx_v, idx_w, idx_p = self.zero_blocks(res, eigs)

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor") and field_col == field_row:
            mat = geo.i4j2(res[0], res[2], m, bc)
            mat = utils.qid_from_idx(idx_u, res[0]*res[2])*mat

        elif field_row == ("velocity","pol") and field_col == field_row:
            mat = geo.i6j4(res[0], res[2], m, bc)
            mat = utils.qid_from_idx(idx_v, res[0]*res[2])*mat

        elif field_row == ("temperature","") and field_col == field_row:
            mat = geo.i2j2(res[0], res[2], m, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        Ta = eq_params['taylor']
        G = eq_params['gamma']
        T = Ta**0.5
        m = eigs[0]

        zscale = 2.0

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                mat = geo.i4j2lapl2(res[0], res[2], m, bc, zscale = zscale)

            elif field_col == ("velocity","pol"):
                mat = geo.i4laplhj2e1(res[0], res[2], m, bc, -T, zscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], res[2], m, 2, 2, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = geo.i6laplhj4e1(res[0], res[2], m, bc, T, zscale = zscale)

            elif field_col == ("velocity","pol"):
                mat = geo.i6j4lapl3(res[0], res[2], m, bc, zscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.i6laplhj4(res[0], res[2], m, bc, -Ra)

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = geo.zblk(res[0], res[2], m, 1, 2, bc)

            elif field_col == ("velocity","pol"):
                if self.linearize or bcs["bcType"] == self.FIELD_TO_RHS:
                    mat = geo.i2laplhj2(res[0], res[2], m, bc, -1.0)
                else:
                    mat = geo.zblk(res[0], res[2], m, 1, 2, bc)

            elif field_col == ("temperature",""):
                mat = geo.i2j2lapl(res[0], res[2], m, bc, 1.0, zscale = zscale)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']
        m = eigs[0]

        zscale = 2.0

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = geo.i4laplhj2(res[0], res[2], m, bc, 1.0/Pr)

        elif field_row == ("velocity","pol"):
            mat = geo.i6j4lapl2(res[0], res[2], m, bc, 1.0/Pr, zscale = zscale)

        elif field_row == ("temperature",""):
            mat = geo.i2j2(res[0], res[2], m, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat
