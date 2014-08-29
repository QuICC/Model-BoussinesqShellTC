"""Module provides the functions to generate the Boussinesq Rayleigh-Benard convection in a 2D box (1 periodic direction) (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_2d import no_bc


class BoussinesqRB2DBoxVC(base_model.BaseModel):
    """Class to setup the Boussinesq Rayleigh-Benard convection in a 2D box (1 periodic direction) (velocity-continuity formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "zxratio"]

    def periodicity(self):
        """Get the domain periodicity"""

        return [False, True, False]

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
            if field_row == ("velocityx","") or field_row == ("velocityy","")  or field_row == ("velocityz","") or field_row == ("pressure","") or field_row == ("temperature",""):
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

    version = 1

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
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

            bc = no_bc()
            bcId = bcs.get(field_col[0], -1)
            # No-slip/No-slip, Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'x'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'x'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'z'}

            # Stress-free/Stress-free, Fixed flux/Fixed flux
            elif bcId == 1:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:21}, 'priority':'x'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:21}, 'z':{0:21}, 'priority':'sx'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'z':{0:21}, 'priority':'sx'}

            # Stress-free/No-slip, Fixed flux/Fixed temperature
            elif bcId == 2:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'x'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'sx'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'z':{0:20}, 'priority':'z'}

            # No-slip/Stress-free, Fixed temperature/Fixed flux
            elif bcId == 3:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'z':{0:21}, 'priority':'x'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:20}, 'z':{0:21}, 'priority':'x'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:20}, 'z':{0:20}, 'priority':'z'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'z':{0:21}, 'priority':'x'}
            
            # Set LHS galerkin restriction
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

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 1:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}

                elif bcId == 2:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 3:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}

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
        return c2d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row):
        """Create the quasi-inverse operator"""

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        elif field_row == ("velocityy",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        elif field_row == ("velocityz",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        elif field_row == ("temperature",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col):
        """Create matrix block linear operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']

        zscale = eq_params['zxratio']

        k = eigs[0]/2.0

        # U: TiN
        zero_ru = spsp.kron(c1d.qid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        # U: Ti0
        zero_cu = spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        idx_ru = (np.ravel(zero_ru.sum(axis=1)) > 0)
        idx_cu = (np.ravel(zero_cu.sum(axis=1)) > 0)
        zero_ru = spsp.lil_matrix(zero_ru.shape)
        zero_cu = spsp.lil_matrix(zero_cu.shape)
        zero_ru[idx_ru,idx_ru] = 1
        zero_cu[idx_cu,idx_cu] = 1
        # U: Handle overlaps T_00 and T_N0
        zero_rut = spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        zero_rut = zero_rut + spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        zero_cu = zero_cu - zero_rut
        idx_rut = (np.ravel(zero_rut.sum(axis=1)) > 0)
        zero_rut = spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-2, c1d.c1dbc.no_bc()))
        idx_cut = (np.ravel(zero_rut.sum(axis=1)) > 0)
        zero_rut = spsp.lil_matrix(zero_rut.shape)
        if self.version == 1:
            zero_rut[idx_cut,idx_rut] = 1
        elif self.version == 2:
            zero_rut[idx_rut,idx_rut] = 1

        # W: TNk
        zero_rw = spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        # W: T0k
        zero_cw = spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        idx_rw = (np.ravel(zero_rw.sum(axis=1)) > 0)
        idx_cw = (np.ravel(zero_cw.sum(axis=1)) > 0)
        zero_rw = spsp.lil_matrix(zero_rw.shape)
        zero_cw = spsp.lil_matrix(zero_cw.shape)
        zero_rw[idx_rw,idx_rw] = 1
        zero_cw[idx_cw,idx_cw] = 1
        # W: Handle overlaps T_00 and T_0N
        zero_rwt = spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        zero_rwt = zero_rwt + spsp.kron(c1d.qid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        zero_cw = zero_cw - zero_rwt
        idx_rwt = (np.ravel(zero_rwt.sum(axis=1)) > 0)
        zero_rwt = spsp.kron(c1d.sid(res[2], res[2]-2, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        idx_cwt = (np.ravel(zero_rwt.sum(axis=1)) > 0)
        zero_rwt = spsp.lil_matrix(zero_rwt.shape)
        if self.version == 1:
            zero_rwt[idx_cwt,idx_rwt] = 1
        elif self.version == 2:
            zero_rwt[idx_rwt,idx_rwt] = 1

        # Pressure: T_iN, T_Nk
        if self.version == 1:
            zero_p = zero_ru + zero_rw
        elif self.version == 2:
            zero_p = zero_ru + zero_rw + zero_cu + zero_cw
        # Pressure: T_{N-2:N,N-2:N}
        zero_p = zero_p + spsp.kron(c1d.qid(res[2], res[2]-3, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-3, c1d.c1dbc.no_bc()))
        # Pressure: T_00
        zero_p = zero_p + spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        idx_p = (np.ravel(zero_p.sum(axis=1)) > 0)
        zero_p = spsp.lil_matrix(zero_p.shape)
        zero_p[idx_p,idx_p] = 1

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = c2d.i2j2lapl(res[0], res[2], k, bc, zscale = zscale).tolil()
                mat[idx_ru,:] = 0;
                if self.version == 2:
                    mat[idx_cu,:] = 0;
                mat[:,idx_ru] = 0;
                mat[:,idx_cu] = 0;
                if self.version == 1:
                    mat = mat + zero_ru + zero_rut
                elif self.version == 2:
                    mat = mat + zero_ru + zero_rut + zero_cu

            elif field_col == ("velocityy",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = c2d.i2j2d1(res[0], res[2], bc, -1.0).tolil()
                mat[idx_ru,:] = 0
                if self.version == 2:
                    mat[idx_cu,:] = 0;
                mat[:,idx_p] = 0

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = c2d.i2j2lapl(res[0], res[2], k, bc, zscale = zscale)

            elif field_col == ("velocityz",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = c2d.i2j2(res[0], res[2], bc, -1j*k).tolil()
                mat[:,idx_p] = 0

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c2d.i2j2lapl(res[0], res[2], k, bc, zscale = zscale).tolil()
                mat[idx_rw,:] = 0;
                if self.version == 2:
                    mat[idx_cw,:] = 0;
                mat[:,idx_rw] = 0;
                mat[:,idx_cw] = 0;
                if self.version == 1:
                    mat = mat + zero_rw + zero_rwt
                elif self.version == 2:
                    mat = mat + zero_rw + zero_rwt + zero_cw

            elif field_col == ("temperature",""):
                mat = c2d.i2j2(res[0], res[2], bc, Ra/16.0).tolil()
                mat[idx_rw,:] = 0;
                if self.version == 2:
                    mat[idx_cw,:] = 0;

            elif field_col == ("pressure",""):
                mat = c2d.i2j2e1(res[0], res[2], bc, -1.0, zscale = zscale).tolil()
                mat[:,idx_p] = 0
                mat[idx_rw,:] = 0
                if self.version == 2:
                    mat[idx_cw,:] = 0;

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c2d.i2j2(res[0], res[2], bc).tolil()
                mat[:,idx_rw] = 0;
                mat[:,idx_cw] = 0;

            elif field_col == ("temperature",""):
                mat = c2d.i2j2lapl(res[0], res[2], 0, bc, zscale = zscale)

            elif field_col == ("pressure",""):
                mat = c2d.zblk(res[0], res[2], 2, 2, bc)

        elif field_row == ("pressure",""):
            if field_col == ("velocityx",""):
                mat = c2d.i1j1d1(res[0], res[2], bc).tolil()
                mat[:,idx_cu] = 0
                mat[idx_cu,:] = 0
                mat[idx_cw,:] = 0
                mat[idx_p,:] = 0
                if self.version == 1:
                    mat = mat + zero_cu

            elif field_col == ("velocityy",""):
                mat = c2d.i1j1(res[0], res[2], bc, 1j*k).tolil()
                mat[idx_cu,:] = 0
                mat[idx_cw,:] = 0
                mat[idx_p,:] = 0

            elif field_col == ("velocityz",""):
                mat = c2d.i1j1e1(res[0], res[2], bc, zscale = zscale).tolil()
                mat[:,idx_cw] = 0
                mat[idx_cu,:] = 0
                mat[idx_cw,:] = 0
                mat[idx_p,:] = 0
                if self.version == 1:
                    mat = mat + zero_cw

            elif field_col == ("temperature",""):
                mat = c2d.zblk(res[0], res[2], 1, 1, bc)

            elif field_col == ("pressure",""):
                mat = c2d.zblk(res[0], res[2], 1, 1, bc).tolil()
                mat = mat + zero_p

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']
        Ra = eq_params['rayleigh']

        k = eigs[0]/2.0

        zero_ru = spsp.kron(c1d.qid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        idx_ru = (np.ravel(zero_ru.sum(axis=1)) > 0)
        zero_cu = spsp.kron(c1d.sid(res[2], res[2]-1, c1d.c1dbc.no_bc()), c1d.qid(res[0], 0, c1d.c1dbc.no_bc()))
        idx_cu = (np.ravel(zero_cu.sum(axis=1)) > 0)
        zero_rw = spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.qid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        idx_rw = (np.ravel(zero_rw.sum(axis=1)) > 0)
        zero_cw = spsp.kron(c1d.qid(res[2], 0, c1d.c1dbc.no_bc()), c1d.sid(res[0], res[0]-1, c1d.c1dbc.no_bc()))
        idx_cw = (np.ravel(zero_cw.sum(axis=1)) > 0)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c2d.i2j2(res[0], res[2], bc, 1.0/Pr).tolil()
            mat[idx_ru,:] = 0
            if self.version == 2:
                mat[idx_cu,:] = 0

        elif field_row == ("velocityy",""):
            mat = c2d.i2j2(res[0], res[2], bc, 1.0/Pr)

        elif field_row == ("velocityz",""):
            mat = c2d.i2j2(res[0], res[2], bc, 1.0/Pr).tolil()
            mat[idx_rw,:] = 0
            if self.version == 2:
                mat[idx_cw,:] = 0

        elif field_row == ("temperature",""):
            mat = c2d.i2j2(res[0], res[2], bc)

        elif field_row == ("pressure",""):
            mat = c2d.zblk(res[0], res[2], 1, 1, bc)

        return mat
