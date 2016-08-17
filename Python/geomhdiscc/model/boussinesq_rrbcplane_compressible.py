"""Module provides the functions to generate the Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Toroidal/poloidal formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as geo
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_1d import no_bc


class BoussinesqRRBCPlaneConfig:
    """Class to setup the Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Toroidal/poloidal formulation)"""

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, True]

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "ekman", "superadiabatic", "adabatic", "polytropic", "gamma", "scale1d", "fast_mean"]

    def automatic_parameters(self, eq_params):
        """Extend parameters with automatically computable values"""

        # Rescale Z direction with ekman number
        d = {"scale1d":eq_params["scale1d"]*eq_params["ekman"]**(1./3.)}

        return d

    def config_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocity", "entropy", "density"]

    def stencil(self, res, eq_params, eigs, bcs, field_row, make_square):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return geo.stencil(res[0], bc, make_square)

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Index mode: SLOWEST_SINGLE_RHS, SLOWEST_MULTI_RHS, MODE, SINGLE
        index_mode = self.MODE

        return self.compile_equation_info(res, field_row, is_complex, index_mode)

class BoussinesqRRBCPlane(BoussinesqRRBCPlaneConfig, base_model.BaseModel):
    """Class to setup the Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Toroidal/poloidal formulation)"""

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocity","x"), ("velocity","y"), ("velocity","z"), ("temperature",""), ("entropy",""), ("density",""), ("pressure","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocity","x"), ("velocity","y"), ("velocity","z"), ("temperature",""), ("entropy",""), ("density",""), ("pressure","")]

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row in [("velocity","x"), ("velocity","y"), ("velocity","z"), ("entropy",""), ("density","")]:
                shift_z = 2
            else:
                shift_z = 0

            gal_n = (res[0] - shift_z)

        else:
            gal_n = tau_n
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_z,0,0), 1)
        return block_info

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            m = eigs[0]

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip / Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocity","x"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("velocity","y"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("entropy",""):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("density",""):
                        bc = {0:-20, 'rt':0}

                else:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("velocity","y") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("entropy","") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("density","") and field_col == field_row:
                        bc = {0:20}

            # Stress-free / Fixed flux
            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocity","x"):
                        bc = {0:-21, 'rt':0}
                    elif field_col == ("velocity","y"):
                        bc = {0:-21, 'rt':0}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20, 'rt':0}
                    elif field_col == ("entropy",""):
                        bc = {0:-21, 'rt':0}
                    elif field_col == ("density",""):
                        bc = {0:-21, 'rt':0}

                else:
                    if field_row == ("velocity","x") and field_col == field_row:
                        bc = {0:21}
                    elif field_row == ("velocity","y") and field_col == field_row:
                        bc = {0:21}
                    elif field_row == ("velocity","z") and field_col == field_row:
                        bc = {0:20}
                    elif field_row == ("entropy","") and field_col == field_row:
                        bc = {0:21}
                    elif field_row == ("density","") and field_col == field_row:
                        bc = {0:21}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row in [("velocity","x"), ("velocity","y"), ("velocity","z")]:
                    bc['rt'] = 2
                elif field_row == ("entropy",""):
                    bc['rt'] = 2
                elif field_row == ("density",""):
                    bc['rt'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocity","x"):
                        bc = {0:-20}
                    elif field_col == ("velocity","y"):
                        bc = {0:-20}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20}
                    elif field_col == ("entropy",""):
                        bc = {0:-20}
                    elif field_col == ("density",""):
                        bc = {0:-20}

                elif bcId == 1:
                    if field_col == ("velocity","x"):
                        bc = {0:-21}
                    elif field_col == ("velocity","y"):
                        bc = {0:-21}
                    elif field_col == ("velocity","z"):
                        bc = {0:-20}
                    elif field_col == ("entropy",""):
                        bc = {0:-21}
                    elif field_col == ("density",""):
                        bc = {0:-21}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row in [("velocity","x"), ("velocity","y"), ("velocity","z")]:
                    bc['rt'] = 2
                elif field_row == ("entropy",""):
                    bc['rt'] = 2
                elif field_row == ("density",""):
                    bc['rt'] = 2

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        kx = eigs[0]
        ky = eigs[1]

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def nonlinear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create the explicit nonlinear operator"""

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def implicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Ra = eq_params['rayleigh']
        Pr = eq_params['prandtl']
        E = eq_params['ekman']
        Ta = E**(-2)
        Ro = E**(1./3.)
        Hs = eq_params['superadiabatic']
        Ha = eq_params['adiabatic']
        poly = eq_params['polytropic']
        gamma = eq_params['gamma']
        zscale = eq_params['scale1d']

        kx = eigs[0]
        ky = eigs[1]

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        # X velocity
        if field_row == ("velocity","x"):
            if field_col == ("velocity","x"):
                mat = geo.i2lapl(res[0], bc, sqrt(Pr/Ra), cscale = zscale)
                bc[0] = min(bc[0], 0)
                mat += geo.i2(res[0], bc, -(1./3.)*sqrt(Pr/Ra)*kx**2)

            elif field_col == ("velocity","y"):
                mat = geo.i2(res[0], bc, sqrt(Pr*Ta/Ra))
                bc[0] = min(bc[0], 0)
                mat += geo.i2(res[0], bc, -(1./3.)*sqrt(Pr/Ra)*kx*ky)

            elif field_col == ("velocity","z"):
                mat += geo.i2d1(res[0], bc, (1./3.)*sqrt(Pr/Ra)*kx*1j, cscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = geo.i2(res[0], bc, -Hs*1j*kx)

        # Y velocity
        elif field_row == ("velocity","y"):
            if field_col == ("velocity","x"):
                mat = geo.i2(res[0], bc, -sqrt(Pr*Ta/Ra))
                bc[0] = min(bc[0], 0)
                mat += geo.i2(res[0], bc, -(1./3.)*sqrt(Pr/Ra)*kx*ky)

            elif field_col == ("velocity","y"):
                mat = geo.i2lapl(res[0], bc, sqrt(Pr/Ra), cscale = zscale)
                bc[0] = min(bc[0], 0)
                mat += geo.i2(res[0], bc, -(1./3.)*sqrt(Pr/Ra)*ky**2)

            elif field_col == ("velocity","z"):
                mat += geo.i2d1(res[0], bc, (1./3.)*sqrt(Pr/Ra)*ky*1j, cscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = geo.i2(res[0], bc, -Hs*ky*1j)

        # Z velocity
        elif field_row == ("velocity","z"):
            if field_col == ("velocity","x"):
                mat = geo.i2d1(res[0], bc, (1./3.)*sqrt(Pr/Ra)*kx*1j, cscale = zscale)

            elif field_col == ("velocity","y"):
                mat += geo.i2d1(res[0], bc, (1./3.)*sqrt(Pr/Ra)*ky*1j, cscale = zscale)

            elif field_col == ("velocity","z"):
                mat = geo.i2lapl(res[0], bc, sqrt(Pr/Ra), cscale = zscale)
                bc[0] = min(bc[0], 0)
                mat += geo.i2d2(res[0], bc, (1./3.)*sqrt(Pr/Ra), cscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.i2(res[0], bc, Hs)

            elif field_col == ("pressure",""):
                mat = geo.i2d1(res[0], bc, -Hs)

        # Density
        elif field_row == ("density",""):
            if field_col == ("velocity","x"):
                mat = geo.i1(res[0], bc, -kx*1j)

            elif field_col == ("velocity","y"):
                mat = geo.i1(res[0], bc, -ky*1j)

            elif field_col == ("velocity","z"):
                mat = geo.i1d1(res[0], bc, -1.0, cscale = zscale)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = geo.zblk(res[0], bc)

        # Entropy
        elif field_row == ("entropy",""):
            if field_col == ("velocity","x"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","y"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","z"):
                mat = geo.i2(res[0], bc, -1.0)

            elif field_col == ("temperature",""):
                mat = geo.i2lapl(res[0], bc, 1.0/sqrt(Pr*Ra), cscale = zscale)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = geo.zblk(res[0], bc)

        # Pressure
        elif field_row == ("pressure",""):
            if field_col == ("velocity","x"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","y"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","z"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = geo.qid(res[0], 0, bc, -1.0)

            elif field_col == ("entropy",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("density",""):
                mat = geo.qid(res[0], 0, bc, -1.0)

            elif field_col == ("pressure",""):
                mat = geo.qid(res[0], 0, bc)

        # Entropy
        elif field_row == ("temperature",""):
            if field_col == ("velocity","x"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","y"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("velocity","z"):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = geo.zblk(res[0], bc)

            elif field_col == ("entropy",""):
                mat = geo.qid(res[0], 0, bc)

            elif field_col == ("density",""):
                mat = geo.qid(res[0], 0, bc)

            elif field_col == ("pressure",""):
                mat = geo.qid(res[0], 0, bc, -1.0/gamma)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        kx = eigs[0]
        ky = eigs[1]

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        # X velocity
        if field_row == ("velocity","x"):
            mat = geo.i2(res[0], bc, 1.0)

        # Y velocity
        elif field_row == ("velocity","y"):
            mat = geo.i2(res[0], bc, 1.0)

        # Z velocity
        elif field_row == ("velocity","z"):
            mat = geo.i2(res[0], bc, 1.0)

        # Entropy
        elif field_row == ("entropy",""):
            mat = geo.i2(res[0], bc, 1.0)

        # Density
        elif field_row == ("density",""):
            mat = geo.i2(res[0], bc, 1.0)

        # Pressure
        elif field_row == ("pressure",""):
            mat = geo.zblk(res[0], bc)

        # Pressure
        elif field_row == ("temperature",""):
            mat = geo.zblk(res[0], bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

    def boundary_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block of boundary operator"""

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        mat = geo.zblk(res[0], bc)
        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat

class BoussinesqRRBCPlaneVisu(BoussinesqRRBCPlaneConfig, base_model.BaseModel):
    """Class to setup the visualization options for TFF scheme """

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  []

        return fields

    def explicit_fields(self, timing, field_row):
        """Get the list of fields with explicit dependence"""

        # Explicit linear terms
        if timing == self.EXPLICIT_LINEAR:
            if field_row in [("mean_temperature", ""),("fluct_temperature", "")]:
                fields = [("temperature","")]
            else:
                fields = []

        # Explicit nonlinear terms
        elif timing == self.EXPLICIT_NONLINEAR:
            fields = []

        # Explicit update terms for next step
        elif timing == self.EXPLICIT_NEXTSTEP:
            fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row in [("velocity","tor"), ("temperature","")]:
                shift_z = 2
            elif field_row in [("velocity","pol")]:
                shift_z = 4
            else:
                shift_z = 0

            gal_n = (res[0] - shift_z)

        else:
            gal_n = tau_n
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_z,0,0), 1)
        return block_info

    def convert_bc(self, eq_params, eigs, bcs, field_row, field_col):
        """Convert simulation input boundary conditions to ID"""

        # Solver: no tau boundary conditions
        if bcs["bcType"] == self.SOLVER_NO_TAU and not self.use_galerkin:
            bc = no_bc()

        # Solver: tau and Galerkin
        elif bcs["bcType"] == self.SOLVER_HAS_BC or bcs["bcType"] == self.SOLVER_NO_TAU:
            raise RuntimeError("Equations are not setup properly!")

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            raise RuntimeError("Equations are not setup properly!")

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                raise RuntimeError("Equations are not setup properly!")

        return bc

    def explicit_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block for explicit linear term"""

        kx = eigs[0]
        ky = eigs[1]

        mat = None
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        # Mean temperature
        if field_row == ("mean_temperature","") and field_col == ("temperature",""):
            if kx == 0 and ky == 0:
                mat = geo.qid(res[0], 0, bc)
            else:
                mat = geo.zblk(res[0], bc)
        elif field_row == ("fluct_temperature","") and field_col == ("temperature",""):
            if kx == 0 and ky == 0:
                mat = geo.zblk(res[0], bc)
            else:
                mat = geo.qid(res[0], 0, bc)

        if mat is None:
            raise RuntimeError("Equations are not setup properly!")

        return mat
