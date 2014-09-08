"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a 1D box (2 periodic directions) (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_1d import no_bc


class BoussinesqRB1DBoxVC(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a 1D box (2 periodic directions) (velocity-continuity formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "zscale"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, True]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocityx", "velocityy", "velocityz", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocityx",""), ("velocityy",""), ("velocityz",""), ("temperature",""), ("pressure","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocityx",""), ("velocityy",""), ("velocityz",""), ("temperature",""), ("pressure","")]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]
        if self.use_galerkin:
            if field_row == ("velocityx","") or field_row == ("velocityy","") or field_row == ("velocityz","") or field_row == ("temperature","") or field_row == ("pressure",""):
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

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: SLOWEST = 0, MODE = 1
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
            m = eigs[0]

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip / Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {0:-20, 'r':0}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {0:-20, 'r':0}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {0:-20, 'r':0}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {0:-20, 'r':0}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {0:20}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {0:20}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {0:20}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {0:20}

            # Stress-free / Fixed flux
            elif bcId == 1:
                if self.use_galerkin:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {0:-21, 'r':0}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {0:-21, 'r':0}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {0:-20, 'r':0}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {0:-21, 'r':0}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {0:21}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {0:21}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {0:20}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {0:21}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_ror == ("velocityz",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityx",""):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("velocityy",""):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("velocityz",""):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-20, 'x':0}

                elif bcId == 1:
                    if field_col == ("velocityx",""):
                        bc = {0:-21, 'x':0}
                    elif field_col == ("velocityy",""):
                        bc = {0:-21, 'x':0}
                    elif field_col == ("velocityz",""):
                        bc = {0:-20, 'x':0}
                    elif field_col == ("temperature",""):
                        bc = {0:-21, 'x':0}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return c1d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c1d.i2(res[0], bc).tolil()
            mat[idx_u,:] = 0;

        elif field_row == ("velocityy",""):
            mat = c1d.i2(res[0], bc).tolil()
            mat[idx_v,:] = 0;

        elif field_row == ("velocityz",""):
            mat = c1d.i2(res[0], bc).tolil()
            mat[idx_w,:] = 0;

        elif field_row == ("temperature",""):
            mat = c1d.i2(res[0], bc)

        elif field_row == ("pressure",""):
            mat = c1d.zblk(res[0], bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']
        zscale = eq_params['zscale']

        k1 = eigs[0]
        k2 = eigs[1]

        zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = c1d.i2lapl(res[0], k1, k2, bc, cscale = zscale).tolil()
                mat[:,idx_u] = 0;
                mat[idx_u,:] = 0
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + zero_u

            elif field_col == ("velocityy",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityz",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = c1d.i2(res[0], bc, -1j*k1).tolil()
                mat[idx_u,:] = 0
                mat[:,idx_p] = 0

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityy",""):
                mat = c1d.i2lapl(res[0], k1, k2, bc, cscale = zscale).tolil()
                mat[:,idx_v] = 0
                mat[idx_v,:] = 0
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + zero_v

            elif field_col == ("velocityz",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("temperature",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("pressure",""):
                mat = c1d.i2(res[0], bc, -1j*k2).tolil()
                mat[:,idx_p] = 0
                mat[idx_v,:] = 0

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityy",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityz",""):
                mat = c1d.i2lapl(res[0], k1, k2, bc, cscale = zscale).tolil()
                mat[:,idx_w] = 0;
                mat[idx_w,:] = 0
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + zero_w 

            elif field_col == ("temperature",""):
                mat = c1d.i2(res[0], bc, Ra)

            elif field_col == ("pressure",""):
                mat = c1d.i2d1(res[0], bc, -zscale).tolil()
                mat[:,idx_p] = 0
                mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityy",""):
                mat = c1d.zblk(res[0], bc)

            elif field_col == ("velocityz",""):
                mat = c1d.i2(res[0], bc).tolil()
                mat[:,idx_w] = 0;

            elif field_col == ("temperature",""):
                mat = c1d.i2lapl(res[0], k1, k2, bc, cscale = zscale)

            elif field_col == ("pressure",""):
                mat = c1d.zblk(res[0], bc)

        elif field_row == ("pressure",""):
            if bcs["bcType"] == self.SOLVER_NO_TAU:
                mat = c1d.zblk(res[0], no_bc())

            else:
                if field_col == ("velocityx",""):
                    bc['rt'] = 1
                    bc['cr'] = 1
                    mat = c1d.i1(res[0]+1, bc, 1j*k1).tolil()
                    mat[:,idx_u] = 0
                    mat[idx_p,:] = 0

                elif field_col == ("velocityy",""):
                    bc['rt'] = 1
                    bc['cr'] = 1
                    mat = c1d.i1(res[0]+1, bc, 1j*k2).tolil()
                    mat[:,idx_v] = 0
                    mat[idx_p,:] = 0

                elif field_col == ("velocityz",""):
                    bc['rt'] = 1
                    bc['cr'] = 1
                    mat = c1d.i1d1(res[0]+1, bc, zscale).tolil()
                    mat[:,idx_w] = 0
                    mat[idx_p,:] = 0

                elif field_col == ("temperature",""):
                    mat = c1d.zblk(res[0], bc)

                elif field_col == ("pressure",""):
                    mat = c1d.zblk(res[0], bc)
                    mat = mat + zero_p

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']

        zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c1d.i2(res[0], bc, 1.0/Pr).tolil()
            mat[:,idx_u] = 0
            mat[idx_u,:] = 0

        elif field_row == ("velocityy",""):
            mat = c1d.i2(res[0], bc, 1.0/Pr).tolil()
            mat[:,idx_v] = 0
            mat[idx_v,:] = 0

        elif field_row == ("velocityz",""):
            mat = c1d.i2(res[0], bc, 1.0/Pr).tolil()
            mat[:,idx_w] = 0
            mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            mat = c1d.i2(res[0], bc)

        elif field_row == ("pressure",""):
            mat = c1d.zblk(res[0], bc)

        return mat

    def zero_blocks(self, res, eigs):
        """Build restriction matrices"""

        # U: TN
        zero_u = c1d.zblk(res[0], no_bc())
        zero_u = zero_u + c1d.qid(res[0], res[0]-1, no_bc())
        # Cleanup and create indexes list
        idx_u = (np.ravel(zero_u.sum(axis=1)) > 0)
        zero_u = spsp.lil_matrix(zero_u.shape)
        zero_u[idx_u,idx_u] = 1

        # V: TN
        zero_v = c1d.zblk(res[0], no_bc())
        zero_v = zero_v + c1d.qid(res[0], res[0]-1, no_bc())
        # Cleanup and create indexes list
        idx_v = (np.ravel(zero_v.sum(axis=1)) > 0)
        zero_v = spsp.lil_matrix(zero_v.shape)
        zero_v[idx_v,idx_v] = 1

        # W:
        zero_w = c1d.zblk(res[0], no_bc())
        # Cleanup and create indexes list
        idx_w = (np.ravel(zero_w.sum(axis=1)) > 0)
        zero_w = spsp.lil_matrix(zero_w.shape)
        zero_w[idx_w,idx_w] = 1

        # Pressure: T_iN, T_Nk
        zero_p = c1d.zblk(res[0], no_bc())
        zero_p = zero_p + c1d.qid(res[0], res[0]-1, no_bc())
        # Cleanup and create indexes list
        idx_p = (np.ravel(zero_p.sum(axis=1)) > 0)
        zero_p = spsp.lil_matrix(zero_p.shape)
        zero_p[idx_p,idx_p] = 1

        return (zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p)
