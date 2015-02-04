"""Module provides the functions to generate the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation) without field coupling (standard implementation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.spherical.shell_radius as shell
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.spherical.shell_radius_boundary import no_bc


class BoussinesqRTCShellStd(base_model.BaseModel):
    """Class to setup the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation) without field coupling (standard implementation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["taylor", "prandtl", "rayleigh", "ro", "rratio", "heating"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  []

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""
    
        #fields = [("velocity","tor"),("velocity","pol"),("temperature","")]
        fields = [field_row]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        if field_row == ("velocity","tor"):
            fields = []
        elif field_row == ("velocity","pol"):
            fields = [("temperature","")]
        elif field_row == ("temperature",""):
            fields = [("velocity","pol")]

        return fields

    def block_size(self, res, field_row):
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

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = False

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.SLOWEST_MULTI_RHS

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

        l = eigs[0]
        a, b = shell.linear_r2x(eq_params['ro'], eq_params['rratio'])

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
                        bc = {0:-40, 'rt':0, 'c':{'a':a, 'b':b}}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':0}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:20}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:40, 'c':{'a':a, 'b':b}}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                            bc = {0:20}

            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':0, 'c':{'a':a, 'b':b}}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':0, 'c':{'a':a, 'b':b}}

                else:
                    if field_row == ("velocity","tor") and field_col == ("velocity","tor"):
                            bc = {0:22, 'c':{'a':a, 'b':b}}
                    elif field_row == ("velocity","pol") and field_col == ("velocity","pol"):
                            bc = {0:41, 'c':{'a':a, 'b':b}}
            
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
                        bc = {0:-40, 'rt':4, 'c':{'a':a, 'b':b}}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':2}

                elif bcId == 1:
                    if field_col == ("velocity","tor"):
                        bc = {0:-22, 'rt':2, 'c':{'a':a, 'b':b}}
                    elif field_col == ("velocity","pol"):
                        bc = {0:-41, 'rt':4, 'c':{'a':a, 'b':b}}
        
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

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return shell.stencil(res[0], bc, make_square)

    def qi(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create the quasi-inverse operator"""

        a, b = shell.linear_r2x(eq_params['ro'], eq_params['rratio'])

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("temperature",""):
            if eq_params["heating"] == 0:
                mat = shell.i2x2(res[0], a, b, bc)
            else:
                mat = shell.i2x3(res[0], a, b, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        Ta = eq_params['taylor']
        ro = eq_params['ro']
        rratio = eq_params['rratio']
        T = Ta**0.5

        # Easy switch from nondimensionalistion by R_o (Dormy) and (R_o - R_i) (Christensen)
        # Parameters match as:  Dormy   Christensen 
        #                       Ra      Ra/(R_o*R_i*Ta^0.5)
        #                       Ta      Ta*(1-R_i/R_o)^4
        if ro == 1.0:
            # R_o rescaling
            Ra_eff = Ra
            bg_eff = 1.0
        elif eq_params['heating'] == 0:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra*T/ro)
            bg_eff = 2.0/(ro*(1.0 + rratio))
        elif eq_params['heating'] == 1:
            # (R_o - R_i) rescaling
            Ra_eff = (Ra*T/ro)
            bg_eff = ro**2*rratio

        l = eigs[0]

        a, b = shell.linear_r2x(eq_params['ro'], eq_params['rratio'])

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocity","tor"):
            if field_col == ("velocity","tor"):
                if l == 0 and bcs['bcType'] == self.SOLVER_HAS_BC:
                    bc[0] = min(0,bc[0]//10)
                    mat = shell.qid(res[0], 0, bc)
                else:
                    mat = shell.i2x2lapl(res[0], l, a, b, bc, l*(l+1.0))

            elif field_col == ("velocity","pol"):
                mat = shell.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = shell.zblk(res[0], bc)

        elif field_row == ("velocity","pol"):
            if field_col == ("velocity","tor"):
                mat = shell.zblk(res[0], bc)

            elif field_col == ("velocity","pol"):
                if l == 0 and bcs['bcType'] == self.SOLVER_HAS_BC:
                    bc[0] = min(0,bc[0]//10)
                    mat = shell.qid(res[0], 0, bc)
                else:
                    mat = shell.i4x4lapl2(res[0], l, a, b, bc, l*(l+1.0))

            elif field_col == ("temperature",""):
                if l == 0 and bcs['bcType'] == self.SOLVER_HAS_BC:
                    bc[0] = min(0,bc[0]//10)
                    mat = shell.zblk(res[0], bc)
                else:
                    mat = shell.i4x4(res[0], a, b, bc, -Ra_eff*l*(l+1.0))

        elif field_row == ("temperature",""):
            if field_col == ("velocity","tor"):
                mat = shell.zblk(res[0], bc)

            elif field_col == ("velocity","pol"):
                if self.linearize or bcs["bcType"] == self.FIELD_TO_RHS:
                    if eq_params["heating"] == 0:
                        mat = shell.i2x2(res[0], a, b, bc, bg_eff*l*(l+1.0))
                    else:
                        mat = shell.i2(res[0], a, b, bc, bg_eff*l*(l+1.0))
             
                else:
                    mat = shell.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                if eq_params["heating"] == 0:
                    mat = shell.i2x2lapl(res[0], l, a, b, bc, 1.0/Pr)
                else:
                    mat = shell.i2x3lapl(res[0], l, a, b, bc, 1.0/Pr)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        l = eigs[0]

        a, b = shell.linear_r2x(eq_params['ro'], eq_params['rratio'])

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocity","tor"):
            mat = shell.i2x2(res[0], a, b, bc, l*(l+1.0))

        elif field_row == ("velocity","pol"):
            mat = shell.i4x4lapl(res[0], l, a, b, bc, l*(l+1.0))

        elif field_row == ("temperature",""):
            if eq_params["heating"] == 0:
                mat = shell.i2x2(res[0], a, b, bc)
            else:
                mat = shell.i2x3(res[0], a, b, bc)

        return mat
