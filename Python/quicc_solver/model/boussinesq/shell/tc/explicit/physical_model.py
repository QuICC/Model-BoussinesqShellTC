"""Module provides the functions to generate the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation) without field coupling (standard implementation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils
import quicc.geometry.spherical.shell_radius as geo
import quicc.base.base_model as base_model
from quicc.geometry.spherical.shell_radius_boundary import no_bc


class PhysicalModel(base_model.BaseModel):
    """Class to setup the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation) without field coupling (standard implementation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "rratio", "heating"]

    def automatic_parameters(self, eq_params):
        """Extend parameters with automatically computable values"""

        d = dict()

        rratio = eq_params['rratio']
        # Unit gap width
        if True:
            gap = {
                    "lower_1d":rratio/(1.0 - rratio),
                    "upper_1d":1.0/(1.0 - rratio)
                    }
        # Unit radius
        else:
            gap = {
                    "lower_1d":rratio,
                    "upper_1d":1.0
                    }

        d.update(gap)

        return d

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","pol"), ("temperature","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields = [field_row]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit linear dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            if field_row == ("velocity","pol"):
                fields = [("temperature","")]
            elif field_row == ("temperature",""):
                fields = [("velocity","pol")]
            else:
                fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            if field_row == ("temperature",""):
                fields = [("temperature","")]
            else:
                fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, eigs, bcs, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row == ("velocity","tor") or field_row == ("temperature",""):
                shift_r = 2
            elif field_row == ("velocity","pol"):
                shift_r = 4
            else:
                shift_r = 0

            gal_n = (res[0] - shift_r)

        else:
            gal_n = tau_n
            shift_r = 0

        block_info = (tau_n, gal_n, (shift_r,0,0), 1)
        return block_info

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""

        assert(eigs[0].is_integer())
        l = eigs[0]

        # Get boundary condition
        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        mat = geo.stencil(res[0], bc, make_square)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is real
        is_complex = False

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_MULTI_RHS

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
                    if field_col == ("velocity","tor"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-40, 'rt':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                        bc = {0:20}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                        bc = {0:40}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {0:20}

            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':0}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                        bc = {0:22}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                        bc = {0:41}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 4
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","tor"):
                        bc = {0:-20, 'rt':2}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-40, 'rt':4}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':2}

                elif bcId == 1:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':2}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':4}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocity","tor"):
                    bc['rt'] = 2
                elif field_row == ("velocity","pol"):
                    bc['rt'] = 4
                elif field_row == ("temperature",""):
                    bc['rt'] = 2

        else:
            bc = no_bc()

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        assert(eigs[0].is_integer())
        l = eigs[0]

        Ra_eff, bg_eff = self.nondimensional_factors(eq_params)

        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","pol") and field_col == ("temperature",""):
            mat = geo.i4r4(res[0], ri, ro, bc, Ra_eff)

        elif field_row == ("temperature","") and field_col == ("velocity","pol"):
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, bc, -bg_eff*l*(l+1.0))
            else:
                mat = geo.i2(res[0], ri, ro, bc, -bg_eff*l*(l+1.0))

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nonlinear term"""

        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == field_row:
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, bc)
            else:
                mat = geo.i2r3(res[0], ri, ro, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        assert(eigs[0].is_integer())
        l = eigs[0]

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']

        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor") and field_col == field_row:
            mat = geo.i2r2lapl(res[0], ri, ro, l, bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","pol"):
                mat = geo.i4r4lapl2(res[0], ri, ro, l, bc)

            elif field_col == ("temperature",""):
                if self.linearize:
                    mat = geo.i4r4(res[0], ri, ro, bc, -Ra)

        elif field_row == ("temperature",""):
            if field_col == ("velocity","pol"):
                if self.linearize:
                    if eq_params["heating"] == 0:
                        mat = geo.i2r2(res[0], ri, ro, bc, l*(l+1.0))
                    else:
                        mat = geo.i2(res[0], ri, ro, bc, l*(l+1.0))

            elif field_col == ("temperature",""):
                if eq_params["heating"] == 0:
                    mat = geo.i2r2lapl(res[0], ri, ro, l, bc, 1.0/Pr)
                else:
                    mat = geo.i2r3lapl(res[0], ri, ro, l, bc, 1.0/Pr)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        assert(eigs[0].is_integer())
        l = eigs[0]

        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = geo.i2r2(res[0], ri, ro, bc)

        elif field_row == ("velocity","pol"):
            mat = geo.i4r4lapl(res[0], ri, ro, l, bc)

        elif field_row == ("temperature",""):
            if eq_params["heating"] == 0:
                mat = geo.i2r2(res[0], ri, ro, bc)
            else:
                mat = geo.i2r3(res[0], ri, ro, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        mat = geo.zblk(res[0], ri, ro, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nondimensional_factors(self, eq_params):
        """Compute the effective Rayleigh number and background depending on nondimensionalisation"""

        Ra = eq_params['rayleigh']
        ri, ro = (self.automatic_parameters(eq_params)['lower_1d'], self.automatic_parameters(eq_params)['upper_1d'])
        rratio = eq_params['rratio']

        # Switch from nondimensionalistion by R_o and (R_o - R_i)
        if ro == 1.0:
            # R_o rescaling
            Ra_eff = Ra
            bg_eff = 1.0
        elif eq_params['heating'] == 0:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra/ro)
            bg_eff = 2.0/(ro*(1.0 + rratio))
        elif eq_params['heating'] == 1:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra/ro)
            bg_eff = ro**2*rratio

        return (Ra_eff, bg_eff)
