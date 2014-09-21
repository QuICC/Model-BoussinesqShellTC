"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a square (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_2d import no_bc


class BoussinesqRBSquareVC(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a square (velocity-continuity formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "scale1d", "scale3d"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, False]

    def all_fields(self):
        """Get the list of fields that need a configuration entry"""

        return ["velocityx", "velocityz", "temperature"]

    def stability_fields(self):
        """Get the list of fields needed for linear stability calculations"""

        fields =  [("velocityx",""), ("velocityz",""), ("temperature",""), ("pressure","")]

        return fields

    def implicit_fields(self, field_row):
        """Get the list of coupled fields in solve"""

        fields =  [("velocityx",""), ("velocityz",""), ("temperature",""), ("pressure","")]

        return fields

    def explicit_fields(self, field_row):
        """Get the list of fields with explicit linear dependence"""

        fields = []

        return fields

    def block_size(self, res, field_row):
        """Create block size information"""

        tau_n = res[0]*res[2]
        if self.use_galerkin:
            if field_row == ("velocityx","") or field_row == ("velocityz","") or field_row == ("temperature","") or field_row == ("pressure",""):
                shift_x = 2
                shift_z = 2
            else:
                shift_x = 0
                shift_z = 0

            gal_n = (res[0] - shift_x)*(res[2] - shift_z)

        else:
            gal_n = tau_n
            shift_x = 0
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_x,0,shift_z), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = True

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode: 
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

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip/No-slip, Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'x'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'z'}

            # Stress-free/No-slip, Fixed flux/Fixed temperature
            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'x'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}

            # Stress-free/Stress-free, Fixed flux/Fixed temperature
            elif bcId == 2:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:21}, 'priority':'sx'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}

            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
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
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 1:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 2:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
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
        return c2d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        zero_u, idx_u, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c2d.i2j2(res[0], res[2], bc).tolil()
            mat[idx_u,:] = 0;

        elif field_row == ("velocityz",""):
            mat = c2d.i2j2(res[0], res[2], bc).tolil()
            mat[idx_w,:] = 0;

        elif field_row == ("temperature",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        elif field_row == ("pressure",""):
            mat = c2d.zblk(res[0], res[2], 1, 1, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Ra = eq_params['rayleigh']

        xscale = eq_params['scale1d']
        zscale = eq_params['scale3d']

        zero_u, idx_u, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = c2d.i2j2lapl(res[0], res[2], 0, bc, xscale = xscale, zscale = zscale).tolil()
                mat[:,idx_u] = 0;
                mat[idx_u,:] = 0
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + zero_u

            elif field_col == ("velocityz",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = c2d.i2j2d1(res[0], res[2], bc, -1.0, xscale = xscale).tolil()
                mat[idx_u,:] = 0
                mat[:,idx_p] = 0

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c2d.i2j2lapl(res[0], res[2], 0, bc, xscale = xscale, zscale = zscale).tolil()
                mat[:,idx_w] = 0;
                mat[idx_w,:] = 0
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + zero_w 

            elif field_col == ("temperature",""):
                mat = c2d.i2j2(res[0], res[2], bc, Ra/16.).tolil()
                mat[idx_w,:] = 0

            elif field_col == ("pressure",""):
                mat = c2d.i2j2e1(res[0], res[2], bc, -1.0, zscale = zscale).tolil()
                mat[:,idx_p] = 0
                mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c2d.i2j2(res[0], res[2], bc).tolil()
                mat[:,idx_w] = 0

            elif field_col == ("temperature",""):
                mat = c2d.i2j2lapl(res[0], res[2], 0, bc, xscale = xscale, zscale = zscale)

            elif field_col == ("pressure",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("pressure",""):
            if bcs["bcType"] == self.SOLVER_HAS_BC:
                if field_col == ("velocityx",""):
                    bc['x']['cr'] = 1
                    bc['x']['rt'] = 1
                    bc['x']['zb'] = 1
                    bc['z']['cr'] = 1
                    bc['z']['rt'] = 1
                    bc['z']['zb'] = 1
                    mat = c2d.i1j1d1(res[0]+1, res[2]+1, bc, xscale = xscale).tolil()
                    mat[:,idx_u] = 0
                    mat[idx_p,:] = 0

                elif field_col == ("velocityz",""):
                    bc['x']['cr'] = 1
                    bc['x']['rt'] = 1
                    bc['x']['zb'] = 1
                    bc['z']['cr'] = 1
                    bc['z']['rt'] = 1
                    bc['z']['zb'] = 1
                    mat = c2d.i1j1e1(res[0]+1, res[2]+1, bc, zscale = zscale).tolil()
                    mat[:,idx_w] = 0
                    mat[idx_p,:] = 0

                elif field_col == ("temperature",""):
                    mat = c2d.zblk(res[0], res[2], 1, 1, bc)

                elif field_col == ("pressure",""):
                    mat = c2d.zblk(res[0], res[2], 1, 1, bc).tolil()
                    mat = mat + zero_p
            else:
                mat = c2d.zblk(res[0], res[2], 1, 1, no_bc())


        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']

        zero_u, idx_u, zero_w, idx_w, zero_p, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c2d.i2j2(res[0], res[2], bc, 1.0/Pr).tolil()
            mat[:,idx_u] = 0
            mat[idx_u,:] = 0

        elif field_row == ("velocityz",""):
            mat = c2d.i2j2(res[0], res[2], bc, 1.0/Pr).tolil()
            mat[:,idx_w] = 0
            mat[idx_w,:] = 0

        elif field_row == ("temperature",""):
            mat = c2d.i2j2(res[0], res[2], bc).tolil()

        elif field_row == ("pressure",""):
            mat = c2d.zblk(res[0], res[2], 1, 1, bc)

        return mat

    def zero_blocks(self, res, eigs):
        """Build restriction matrices"""

        # U: TiN
        zero_u = c2d.zblk(res[0], res[2], 2, 2, no_bc())
        zero_u = zero_u + spsp.kron(c1d.qid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        # Cleanup and create indexes list
        idx_u = (np.ravel(zero_u.sum(axis=1)) > 0)
        zero_u = spsp.lil_matrix(zero_u.shape)
        zero_u[idx_u,idx_u] = 1

        # W: TNk
        zero_w = c2d.zblk(res[0], res[2], 2, 2, no_bc())
        zero_w = zero_w + spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        # Cleanup and create indexes list
        idx_w = (np.ravel(zero_w.sum(axis=1)) > 0)
        zero_w = spsp.lil_matrix(zero_w.shape)
        zero_w[idx_w,idx_w] = 1

        # Pressure: T_iN, T_Nk
        zero_p = c2d.zblk(res[0], res[2], 2, 2, no_bc())
        zero_p = zero_p + spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        zero_p = zero_p + spsp.kron(c1d.qid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        # Pressure: T_{N-2:N,N-2:N}
        zero_p = zero_p + spsp.kron(c1d.qid(res[2], res[2]-3, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-3, c1d.c1dbc.no_bc()))
        # Pressure: T_00
        zero_p = zero_p + spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        # Cleanup and create indexes list
        idx_p = (np.ravel(zero_p.sum(axis=1)) > 0)
        zero_p = spsp.lil_matrix(zero_p.shape)
        zero_p[idx_p,idx_p] = 1

        return (zero_u, idx_u, zero_w, idx_w, zero_p, idx_p)