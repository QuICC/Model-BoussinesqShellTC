"""Module provides the functions to generate the periodic Boussinesq Beta 3DQG model"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_1d import no_bc


class BoussinesqBeta3DQGPer(base_model.BaseModel):
    """Class to setup the periodic Boussinesq Beta 3DQG model"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "gamma", "chi", "scale1d"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, True]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["streamfunction", "velocityz", "temperature", "vorticityz"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields = [("streamfunction",""), ("velocityz",""), ("temperature",""), ("vorticityz", "")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        if field_row == ("streamfunction","") or field_row == ("velocityz","") or field_row == ("temperature","") or field_row == ("vorticityz",""):
            fields = [("streamfunction",""), ("velocityz",""), ("temperature",""), ("vorticityz","")]
        else:
            fields = []

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row == ("temperature","") or field_row == ("velocityz",""):
                shift_x = 0
            elif field_row == ("streamfunction",""):
                shift_x = 0
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

        # Matrix operator is complex
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: 
        index_mode = self.MODE

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
            tanchi = np.tan(eq_params['chi']*np.pi/180)
            G = eq_params['gamma']
            kx = eigs[0]
            ky = eigs[1]

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            if bcId == 0:
                if self.use_galerkin:
                    if field_row == ("streamfunction","") and field_col == ("streamfunction",""):
                        if bcs["bcType"] == self.SOLVER_NO_TAU:
                            bc = {0:0}
                        else:
                            bc = {0:10, 'c':-1j*ky*tanchi/G}
                    elif field_row == ("streamfunction","") and field_col == ("velocityz",""):
                        if bcs["bcType"] == self.SOLVER_NO_TAU:
                            bc = {0:0}
                        else:
                            bc = {0:10}
                    elif field_row == ("velocityz","") and field_col == ("streamfunction",""):
                        if bcs["bcType"] == self.SOLVER_NO_TAU:
                            bc = {0:0}
                        else:
                            bc = {0:11, 'c':1j*ky*tanchi/G}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        if bcs["bcType"] == self.SOLVER_NO_TAU:
                            bc = {0:0}
                        else:
                            bc = {0:11}

                else:
                    if field_row == ("vorticityz","") and field_col == ("streamfunction",""):
                        bc = {0:10, 'c':-1j*ky*tanchi/G}
                    elif field_row == ("vorticityz","") and field_col == ("velocityz",""):
                        bc = {0:10}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {0:11}
                    elif field_row == ("velocityz","") and field_col == ("streamfunction",""):
                        bc = {0:11, 'c':1j*ky*tanchi/G}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("streamfunction",""):
                    bc['rt'] = 0
                elif field_row == ("velocityz",""):
                    bc['rt'] = 0
                elif field_row == ("temperature",""):
                    bc['rt'] = 0
                elif field_row == ("vorticity",""):
                    bc['rt'] = 0

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("streamfunction",""):
                        bc = {0:0}
                    elif field_col == ("velocityz",""):
                        bc = {0:0}
                    elif field_col == ("temperature",""):
                        bc = {0:0}
        
        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()

        else:
            bc = no_bc()

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return c1d.stencil(res[0], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create the quasi-inverse operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityz",""):
            mat = c1d.i1(res[0], bc)

        elif field_row == ("vorticityz",""):
            mat = c1d.i1(res[0], bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block of linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        G = eq_params['gamma']
        tanchi = np.tan(eq_params['chi']*np.pi/180)
        zscale = eq_params['scale1d']

        kx = eigs[0]
        ky = eigs[1]

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("streamfunction",""):
            if field_col == ("streamfunction",""):
                mat = c1d.qid(res[0], 0, bc, (kx**2 + ky**2))

            elif field_col == ("velocityz",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("vorticityz",""):
                mat = c1d.qid(res[0], 0, bc)

        elif field_row == ("velocityz",""):
            if field_col == ("streamfunction",""):
                mat = c1d.i1d1(res[0], bc, (-1.0/G**2), cscale = zscale)

            elif field_col == ("velocityz",""):
                mat = c1d.i1(res[0], bc, -(kx**2 + ky**2))

            elif field_col == ("temperature",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("vorticityz",""):
                mat = c1d.zblk(res[0], bc)

        elif field_row == ("temperature",""):
            if field_col == ("streamfunction",""):
                if self.linearize:
                    mat = c1d.qid(res[0], 0, bc, 1j*ky)
                else:
                    mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityz",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = c1d.qid(res[0], 0, bc, -(1/Pr)*(kx**2 + ky**2))

            elif field_col == ("vorticityz",""):
                mat = c1d.zblk(res[0], bc)

        elif field_row == ("vorticityz",""):
            if field_col == ("streamfunction",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityz",""):
                mat = c1d.i1d1(res[0], bc, cscale = zscale)

            elif field_col == ("temperature",""):
                mat = c1d.i1(res[0], bc, 1j*ky*(Ra/(16.0*Pr)))

            elif field_col == ("vorticityz",""):
                mat = c1d.i1(res[0], bc, -(kx**2 + ky**2))

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("streamfunction",""):
            mat = c1d.zblk(res[0], bc)

        elif field_row == ("velocityz",""):
            mat = c1d.i1(res[0], bc)

        elif field_row == ("temperature",""):
            mat = c1d.qid(res[0], 0, bc)

        elif field_row == ("vorticityz",""):
            mat = c1d.i1(res[0], bc)

        return mat
