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

        return ["prandtl", "rayleigh", "rratio", "alpha", "beta", "gamma"]

    def automatic_parameters(self, eq_params):
        """Extend parameters with automatically computable values"""

        d = dict()

        rratio = eq_params['r_ratio']

        useGapWidth = False;
        # Unit gap width
        if useGapWidth:
            gap = {
                    "lower1d":rratio/(1. - rratio),
                    "upper1d":1.0/(1. - rratio)
                    }
        # Unit radius
        else:
            gap = {
                    "lower1d":rratio,
                    "upper1d":1.0
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
        ll1 = l*(l + 1)

        # A = 1 yields a linear gravity in r, whereas A = 0 leads to a 1/r^2 gravity profile
        # B = 1 yields a linear background temperature gradient in r, whereas B = 0 leads to a more general temperature gradient profile
	    # G = 0 solves the linear onset problem; G =/= 0 solves the eigenproblem associated with the energy method
        Ra = eq_params['rayleigh']
        Pr = eq_params['prandtl']
        A = eq_params['alpha']
        B = eq_params['beta']
        G = eq_params['gamma']

        Ra_eff = (Ra/Pr)*ll1

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if G == 0:
            if field_row == ("velocity","pol") and field_col == ("temperature",""):
                mat = geo.i4r4(res[0], ri, ro, bc, Ra_eff*A) + geo.i4r1(res[0], ri, ro, bc, Ra_eff*(1. - A))

            elif field_row == ("temperature","") and field_col == ("velocity","pol"):
                if B == 1:
                    mat = geo.i2r2(res[0], ri, ro, bc, -ll1*B)
                else:
                    mat = geo.i2r3(res[0], ri, ro, bc, -ll1*B) + geo.i2(res[0], ri, ro, bc, -ll1*(1. - B))
        else:
            if field_row == ("velocity","pol") and field_col == ("temperature",""):
                c1 = 0.5*Ra_eff*(A + G**2*B)
                c2 = 0.5*Ra_eff*((1. - A)+G**2*(1. - B))
                mat = geo.i4r4(res[0], ri, ro, bc, c1) + geo.i4r1(res[0], ri, ro, bc, c2) 

            elif field_row == ("temperature","") and field_col == ("velocity","pol"):
                if A == 1 and B == 1:
                    mat = geo.i2r2(res[0], ri, ro, bc, -ll1*0.5*(A/G**2 + B))
                else:
                    c1 = -ll1*0.5*(A/G**2 + B)
                    c2 = -ll1*0.5*((1. - A) + G**2*(1. - B))/G**2
                    mat = geo.i2r3(res[0], ri, ro, bc, c1) + geo.i2(res[0], ri, ro, bc, c2)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit nonlinear term"""

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        # A = 1 yields a linear gravity in r, whereas A = 0 leads to a 1/r^2 gravity profile
        # B = 1 yields a linear background temperature gradient in r, whereas B = 0 leads to a more general temperature gradient profile
	    # G = 0 solves the linear onset problem; G =/= 0 solves the eigenproblem associated with the energy method
        A = eq_params['alpha']
        B = eq_params['beta']
        G = eq_params['gamma']

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("temperature","") and field_col == field_row:
            if G == 0:
                if B == 1:
                    mat = geo.i2r2(res[0], ri, ro, bc)
                else:
                    mat = geo.i2r3(res[0], ri, ro, bc)
             else:
                if A == 1 and B == 1:
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
        ll1 = l*(l + 1.)

        # A = 1 yields a linear gravity in r, whereas A = 0 leads to a 1/r^2 gravity profile
        # B = 1 yields a linear background temperature gradient in r, whereas B = 0 leads to a more general temperature gradient profile
	    # G = 0 solves the linear onset problem; G =/= 0 solves the eigenproblem associated with the energy method
        Ra = eq_params['rayleigh']
        Pr = eq_params['prandtl']
        A = eq_params['alpha']
        B = eq_params['beta']
        G = eq_params['gamma']

        Ra_eff = Ra/Pr*ll1

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if G == 0:
            if field_row == ("velocity","tor") and field_col == field_row:
                mat = geo.i2r2lapl(res[0], ri, ro, l, bc, ll1)

            elif field_row == ("velocity","pol"):
                if field_col == ("velocity","pol"):
                    mat = geo.i4r4lapl2(res[0], ri, ro, l, bc, ll1)

                elif field_col == ("temperature",""):
                    if self.linearize:
                        c1 = -Ra_eff*A
                        c2 = -Ra_eff*(1. - A)
                        mat = geo.i4r4(res[0], ri, ro, bc, c1) + geo.i4r1(res[0], ri, ro, bc, c2)

            elif field_row == ("temperature",""):
                if field_col == ("velocity","pol"):
                    if self.linearize:
                        if B == 1:
                            mat = geo.i2r2(res[0], ri, ro, bc, ll1*B)
                        else:
                            c1 = ll1*B
                            c2 = ll1*(1. - B)
                            mat = geo.i2r3(res[0], ri, ro, bc, c1) + geo.i2(res[0], ri, ro, bc, c2)

                elif field_col == ("temperature",""):
                    if B == 1:
                        mat = geo.i2r2lapl(res[0], ri, ro, l, bc, 1.0/Pr)
                    else:
                        mat = geo.i2r3lapl(res[0], ri, ro, l, bc, 1.0/Pr)
        else:
            if field_row == ("velocity","tor") and field_col == field_row:
                mat = geo.i2r2lapl(res[0], ri, ro, l, bc, ll1)

            elif field_row == ("velocity","pol"):
                if field_col == ("velocity","pol"):
                    mat = geo.i4r4lapl2(res[0], ri, ro, l, bc, ll1)

                elif field_col == ("temperature",""):
                    if self.linearize:
                        c1 = -0.5*Ra_eff*(A + G**2*B)
                        c2 = -0.5*Ra_eff*((1. - A)+G**2*(1 - B))
                        mat = geo.i4r4(res[0], ri, ro, bc, c1) + geo.i4r1(res[0], ri, ro, bc, c2)

            elif field_row == ("temperature",""):
                if field_col == ("velocity","pol"):
                    if self.linearize:
                        if A == 1 and B == 1:
                            mat = geo.i2r2(res[0], ri, ro, bc, 0.5*ll1*(A/G**2 + B))
                        else:
                            c1 = 0.5*ll1*(A/G**2 + B)
                            c2 = 0.5*ll1*((1. - A)+G**2*(1. - B))/G**2
                            mat = geo.i2r3(res[0], a, b, bc, c1) + geo.i2(res[0], ri, ro, bc, c2)

                elif field_col == ("temperature",""):
                    if A == 1 and B == 1:
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
        ll1 = l*(l + 1.)

        # A = 1 yields a linear gravity in r, whereas A = 0 leads to a 1/r^2 gravity profile
        # B = 1 yields a linear background temperature gradient in r, whereas B = 0 leads to a more general temperature gradient profile
	    # G = 0 solves the linear onset problem; G =/= 0 solves the eigenproblem associated with the energy method
        A = eq_params['alpha']
        B = eq_params['beta']
        G = eq_params['gamma']

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if G == 0:
            if field_row == ("velocity","tor"):
                mat = geo.i2r2(res[0], ri, ro, bc, ll1)

            elif field_row == ("velocity","pol"):
                mat = geo.i4r4lapl(res[0], ri, ro, l, bc, ll1)

            elif field_row == ("temperature",""):
                if B == 1:
                    mat = geo.i2r2(res[0], ri, ro, bc)
                else:
                    mat = geo.i2r3(res[0], ri, ro, bc)
        else:
            if field_row == ("velocity","tor"):
                mat = geo.i2r2(res[0], ri, ro, bc, ll1)

            elif field_row == ("velocity","pol"):
                mat = geo.i4r4lapl(res[0], ri, ro, l, bc, ll1)

            elif field_row == ("temperature",""):
                if A == 1 and B == 1:
                    mat = geo.i2r2(res[0], ri, ro, bc)
                else:
                    mat = geo.i2r3(res[0], ri, ro, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        ri, ro = (self.automatic_parameters(eq_params)['lower1d'], self.automatic_parameters(eq_params)['upper1d'])

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        mat = geo.zblk(res[0], ri, ro, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat