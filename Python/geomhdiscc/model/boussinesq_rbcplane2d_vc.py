"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a plane layer (2D) (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as geo
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_1d import no_bc


class BoussinesqRBCPlane2DVC(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a plane layer (2D) (velocity-continuity formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "heating", "scale1d"]

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","x"), ("velocity","z"), ("temperature",""), ("pressure","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocity","x"), ("velocity","z"), ("temperature",""), ("pressure","")]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            if field_row == ("temperature",""):
                fields = [("velocity","z")]
            else:
                fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row in [("velocity","x"), ("velocity","z"), ("temperature","")]:
                fields = [field_row]
            else:
                fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row in [("velocity","x"), ("velocity","z"), ("temperature",""), ("pressure","")]:
                shift_x = 2
            else:
                shift_x = 0

            gal_n = (res[0] - shift_x)

        else:
            gal_n = tau_n
            shift_x = 0

        block_info = (tau_n, gal_n, (shift_x,0,0), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE 
        index_mode = self.MODE

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
            # No-slip / Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:-20, 'rt':0}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:-20, 'rt':0}
                    elif field_row == ("temperature","") and field_col == field_row:
                        bc = {0:-20, 'rt':0}

                else:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("temperature","") and field_col == field_row:
                        bc = {0:20}

            # Stress-free / Fixed flux
            elif bcId == 1:
                if self.use_galerkin:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:-21, 'rt':0}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:-20, 'rt':0}
                    elif field_row == ("temperature","") and field_col == field_row:
                        bc = {0:-21, 'rt':0}

                else:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:21}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("temperature","") and field_col == field_row:
                        bc = {0:21}

            # Fixed temperature at top/Fixed flux at bottom
            elif bcId == 2:
                if self.use_galerkin:
                    if field_row == ("temperature","") and field_col == field_row:
                        bc = {0:-22, 'rt':0}

                else:
                    if field_row == ("temperature","") and field_col == field_row:
                        bc = {0:22}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","x"):
                    bc['rt'] = 2
                elif field_ror == ("velocity","z"):
                    bc['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","x"):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'x':0}

                elif bcId == 1:
                    if field_col == ("velocity","x"):
                        bc = {0:-21, 'x':0}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-21, 'x':0}

                elif bcId == 2:
                    if field_col == ("temperature",""):
                        bc = {0:-22, 'x':0}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","x"):
                    bc['rt'] = 2
                elif field_row == ("velocity","z"):
                    bc['rt'] = 2
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == ("velocity","z"):
            if eq_params['heating'] == 0:
                mat = geo.i2(res[0], bc, -1.0)

            elif eq_params['heating'] == 1:
                mat = geo.i2x1(res[0], bc, -1.0)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","x") and field_col == field_row:
            mat = geo.i2(res[0], bc)

        elif field_row == ("velocity","z") and field_col == field_row:
            mat = geo.i2(res[0], bc)

        elif field_row == ("temperature","") and field_col == field_row:
            mat = geo.i2(res[0], bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Ra = eq_params['rayleigh']
        zscale = eq_params['scale1d']
        
        k = eigs[0]

        idx_p = self.zero_blocks(res, eigs)

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","x"):
            if field_col == ("velocity","x"):
                mat = geo.i2lapl(res[0], k, 0, bc, cscale = zscale)

            elif field_col == ("velocity","z"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = geo.i2(res[0], bc, -1j*k)*utils.qid_from_idx(idx_p, res[0])

        elif field_row == ("velocity","z"):
            if field_col == ("velocity","x"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","z"):
                mat = geo.i2lapl(res[0], k, 0, bc, cscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.i2(res[0], bc, Ra)

            elif field_col == ("pressure",""):
                mat = geo.i2d1(res[0], bc, -1.0, cscale = zscale)*utils.qid_from_idx(idx_p, res[0])

        elif field_row == ("temperature",""):
            if field_col == ("velocity","x"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","z"):
                if self.linearize or bcs["bcType"] == self.FIELD_TO_RHS:
                    if eq_params['heating'] == 0:
                        mat = geo.i2(res[0], bc)

                    elif eq_params['heating'] == 1:
                        mat = geo.i2x1(res[0], bc)
                else:
                    mat = geo.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = geo.i2lapl(res[0], k, 0, bc, cscale = zscale)

            elif field_col == ("pressure",""):
                mat = geo.zblk(res[0], bc)

        elif field_row == ("pressure",""):
            if bcs["bcType"] == self.SOLVER_NO_TAU:
                mat = geo.zblk(res[0], no_bc())

            else:
                if field_col == ("velocity","x"):
                    mat = geo.i2(res[0], bc, 1j*k)

                elif field_col == ("velocity","z"):
                    mat = geo.i2d1(res[0], bc, cscale = zscale)

                elif field_col == ("temperature",""):
                    mat = geo.zblk(res[0], bc)

                elif field_col == ("pressure",""):
                    mat = geo.zblk(res[0], bc)
                    zero_p = spsp.lil_matrix(mat.shape)
                    zero_p[0, -2] = 1
                    zero_p[1, -1] = 1
                    mat = mat + zero_p

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","x"):
            mat = geo.i2(res[0], bc, 1.0/Pr)

        elif field_row == ("velocity","z"):
            mat = geo.i2(res[0], bc, 1.0/Pr)

        elif field_row == ("temperature",""):
            mat = geo.i2(res[0], bc)

        elif field_row == ("pressure",""):
            mat = geo.zblk(res[0], bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def zero_blocks(self, res, eigs):
        """Build restriction matrices"""

        # Pressure: T_N
        idx_p = utils.qidx(res[0], res[0]-2)

        return idx_p
