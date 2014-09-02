"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a cylindrical annulus (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cylindrical.annulus as annulus
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cylindrical.annulus_boundary import no_bc


class BoussinesqRBAnnulusVC(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a cylindrical annulus (velocity-continuity formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "ro", "rratio", "zrratio"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, False, False]

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

        tau_n = res[0]*res[2]
        if self.use_galerkin:
            if field_row == ("velocityx","") or field_row == ("velocityy","") or field_row == ("velocityz","") or field_row == ("temperature",""):
                shift_r = 2
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

        # Matrix operator is complex
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: SLOWEST = 0, MODE = 1
        index_mode = self.SLOWEST

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
            a, b = annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip/No-slip, Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'r':{0:20}, 'z':{0:20}, 'priority':'r'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'r':{0:20}, 'z':{0:20}, 'priority':'r'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'r':{0:20}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'r':{0:20}, 'z':{0:20}, 'priority':'r'}

            # Stress-free/Stress-free, Fixed flux/Fixed flux
            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-24, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'r':{0:20}, 'z':{0:21}, 'priority':'r'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'r':{0:24, 'c':{'a':a, 'b':b}}, 'z':{0:21}, 'priority':'sr'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'r':{0:21}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'r':{0:21}, 'z':{0:21}, 'priority':'sr'}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityx",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 1:
                    if field_col == ("velocityx",""):
                        bc = {'r':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'r':{0:-24, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'r':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'r':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['r']['r'] = 2
                    bc['z']['r'] = 2

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return c2d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        a, b = annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc)

        elif field_row == ("velocityy",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc)

        elif field_row == ("velocityz",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc)

        elif field_row == ("temperature",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']

        zscale = 2.0*eq_params['zrratio']

        a, b = annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])
        m = eigs[0]

        zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = annulus.i2j2x2lapl(res[0], res[2], m, a, b, bc, zscale = zscale)
                bc['r'][0] = min(bc['r'][0], 0)
                bc['z'][0] = min(bc['z'][0], 0)
                mat = mat + annulus.i2j2(res[0], res[2], a, b, bc, -1.0).tolil()
                mat[:,idx_u] = 0
                mat[idx_u,:] = 0

            elif field_col == ("velocityy",""):
                mat = annulus.i2j2(res[0], res[2], a, b, bc, -2.0*1j*m).tolil()
                mat[idx_u,:] = 0

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x2d1(res[0], res[2], a, b, bc, -1.0).tolil()
                mat[idx_u,:] = 0

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = annulus.i2j2(res[0], res[2], a, b, bc, 2.0*1j*m).tolil()
                mat[idx_v,:] = 0

            elif field_col == ("velocityy",""):
                mat = annulus.i2j2x2lapl(res[0], res[2], m, a, b, bc, zscale = zscale)
                bc['r'][0] = min(bc['r'][0], 0)
                bc['z'][0] = min(bc['z'][0], 0)
                mat = mat + annulus.i2j2(res[0], res[2], a, b, bc, -1.0).tolil()
                mat[:,idx_v] = 0
                mat[idx_v,:] = 0

            elif field_col == ("velocityz",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x1(res[0], res[2], a, b, bc, -1j*m).tolil()
                mat[:,idx_p] = 0
                mat[idx_v,:] = 0

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = annulus.i2j2x2lapl(res[0], res[2], m, a, b, bc, zscale = zscale).tolil()
                mat[:,idx_w] = 0
                mat[idx_w,:] = 0

            elif field_col == ("temperature",""):
                mat = annulus.i2j2x2(res[0], res[2], a, b, bc, Ra).tolil()
                mat[idx_w,:] = 0

            elif field_col == ("pressure",""):
                mat = annulus.i2j2x2e1(res[0], res[2], a, b, bc, -1.0, zscale = zscale).tolil()
                mat[:,idx_p] = 0
                mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = annulus.i2j2x2(res[0], res[2], a, b, bc).tolil()
                mat[:,idx_w] = 0

            elif field_col == ("temperature",""):
                mat = annulus.i2j2x2lapl(res[0], res[2], m, a, b, bc, zscale = zscale)

            elif field_col == ("pressure",""):
                mat = annulus.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("pressure",""):
            if field_col == ("velocityx",""):
                mat = annulus.i1j1x1div(res[0], res[2], a, b, bc).tolil()
                mat[idx_p,:] = 0

            elif field_col == ("velocityy",""):
                mat = annulus.i1j1(res[0], res[2], a, b, bc, 1j*m).tolil()
                mat[idx_p,:] = 0

            elif field_col == ("velocityz",""):
                mat = annulus.i1j1x1e1(res[0], res[2], a, b, bc, zscale = zscale).tolil()
                mat[idx_p,:] = 0

            elif field_col == ("temperature",""):
                mat = annulus.zblk(res[0], res[2], 1, 1, bc)

            elif field_col == ("pressure",""):
                mat = annulus.zblk(res[0], res[2], 1, 1, bc)

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']

        a, b = annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])
        m = eigs[0]

        zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc, 1.0/Pr).tolil()
            mat[:,idx_u] = 0
            mat[idx_u,:] = 0

        elif field_row == ("velocityy",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc, 1.0/Pr).tolil()
            mat[:,idx_v] = 0
            mat[idx_v,:] = 0

        elif field_row == ("velocityz",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc, 1.0/Pr).tolil()
            mat[:,idx_w] = 0
            mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            mat = annulus.i2j2x2(res[0], res[2], a, b, bc)

        elif field_row == ("pressure",""):
            mat = annulus.zblk(res[0], res[2], 2, 2, bc)

        return mat

    def zero_blocks(self, res, eigs):
        """Build restriction matrices"""

        # U:
        zero_u = annulus.zblk(res[0], res[1], 2, 2, no_bc())
        # Cleanup and create indexes list
        idx_u = (np.ravel(zero_u.sum(axis=1)) > 0)
        zero_u = spsp.lil_matrix(zero_u.shape)
        zero_u[idx_u,idx_u] = 1

        # V:
        zero_v = annulus.zblk(res[0], res[1], 1, 2, no_bc())
        # Cleanup and create indexes list
        idx_v = (np.ravel(zero_v.sum(axis=1)) > 0)
        zero_v = spsp.lil_matrix(zero_v.shape)
        zero_v[idx_v,idx_v] = 1

        # W:
        zero_w = annulus.zblk(res[0], res[1], 2, 2, no_bc())
        # Cleanup and create indexes list
        idx_w = (np.ravel(zero_w.sum(axis=1)) > 0)
        zero_w = spsp.lil_matrix(zero_w.shape)
        zero_w[idx_w,idx_w] = 1

        # P:
        zero_p = annulus.zblk(res[0], res[1], 2, 2, no_bc())
        # Cleanup and create indexes list
        idx_p = (np.ravel(zero_p.sum(axis=1)) > 0)
        zero_p = spsp.lil_matrix(zero_p.shape)
        zero_p[idx_p,idx_p] = 1

        return (zero_u, idx_u, zero_v, idx_v, zero_w, idx_w, zero_p, idx_p)
