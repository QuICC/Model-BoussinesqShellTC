"""Module provides the functions to generate the Boussinesq rotating Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import geomhdiscc.base.utils as utils
import geomhdiscc.geometry.cartesian.cartesian_3d as c3d
import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
import geomhdiscc.base.base_model as base_model
from geomhdiscc.geometry.cartesian.cartesian_boundary_3d import no_bc


class BoussinesqRRB3DBoxVC(base_model.BaseModel):
    """Class to setup the Boussinesq rotating Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)"""

    def nondimensional_parameters(self):
        """Get the list of nondimensional parameters"""

        return ["prandtl", "rayleigh", "taylor", "scale1d", "scale2d", "scale3d"]

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

        tau_n = res[0]*res[1]*res[2]
        if self.use_galerkin:
            if field_row == ("velocityx","") or field_row == ("velocityy","")  or field_row == ("velocityz","") or field_row == ("temperature","") or field_row == ("pressure",""):
                shift_x = 2
                shift_y = 2
                shift_z = 2
            else:
                shift_x = 0
                shift_y = 0
                shift_z = 0

            gal_n = (res[0] - shift_x)*(res[1] - shift_y)*(res[2] - shift_z)

        else:
            gal_n = tau_n
            shift_x = 0
            shift_y = 0
            shift_z = 0

        block_info = (tau_n, gal_n, (shift_x,shift_y,shift_z), 1)
        return block_info

    def equation_info(self, res, field_row):
        """Provide description of the system of equation"""

        # Matrix operator is complex except for vorticity and mean temperature
        is_complex = False

        # Implicit field coupling
        im_fields = self.implicit_fields(field_row)
        # Additional explicit linear fields
        ex_fields = self.explicit_fields(field_row)

        # Index mode:
        index_mode = self.SINGLE

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
            # No-slip/No-slip/No-slip, Fixed temperature/Fixed temperature/Fixed temperature
            if bcId == 0:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'xz'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'yx'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'zy'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:20}, 'y':{0:20}, 'z':{0:20}, 'priority':'zy'}

            # Stress-free/Stress-free/no-slip, Fixed flux/Fixed flux/Fixed temperature
            elif bcId == 4:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'y':{0:21}, 'z':{0:20}, 'priority':'xz'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:21}, 'y':{0:20}, 'z':{0:20}, 'priority':'yz'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'y':{0:21}, 'z':{0:20}, 'priority':'zsx'}
                    elif field_row == ("temperature","") and field_col == ("temperature",""):
                        bc = {'x':{0:21}, 'y':{0:21}, 'z':{0:20}, 'priority':'zsy'}

            # Stress-free/Stress-free/no-slip, Fixed flux/Fixed flux/Fixed temperature
            elif bcId == 6:
                if self.use_galerkin:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-21, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

                else:
                    if field_row == ("velocityx","") and field_col == ("velocityx",""):
                        bc = {'x':{0:20}, 'y':{0:21}, 'z':{0:21}, 'priority':'xsz'}
                    elif field_row == ("velocityy","") and field_col == ("velocityy",""):
                        bc = {'x':{0:21}, 'y':{0:20}, 'z':{0:21}, 'priority':'ysz'}
                    elif field_row == ("velocityz","") and field_col == ("velocityz",""):
                        bc = {'x':{0:21}, 'y':{0:21}, 'z':{0:20}, 'priority':'zsx'}
            
            # Set LHS galerkin restriction
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2

        # Stencil:
        elif bcs["bcType"] == self.STENCIL:
            if self.use_galerkin:
                bcId = bcs.get(field_col[0], -1)
                if bcId == 0:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}

                elif bcId == 4:
                    if field_col == ("velocityx",""):
                        bc = {'x':{0:-20, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityy",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-20, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("velocityz",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}
                    elif field_col == ("temperature",""):
                        bc = {'x':{0:-21, 'r':0}, 'y':{0:-21, 'r':0}, 'z':{0:-20, 'r':0}}

        # Field values to RHS:
        elif bcs["bcType"] == self.FIELD_TO_RHS:
            bc = no_bc()
            if self.use_galerkin:
                if field_row == ("velocityx",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityy",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("velocityz",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2
                elif field_row == ("temperature",""):
                    bc['x']['r'] = 2
                    bc['y']['r'] = 2
                    bc['z']['r'] = 2

        return bc

    def stencil(self, res, eq_params, eigs, bcs, field_row):
        """Create the galerkin stencil"""
        
        # Get boundary condition
        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        return c3d.stencil(res[0], res[2], bc)

    def qi(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create the quasi-inverse operator"""

        idx_u, idx_v, idx_w, idx_p = self.zero_blocks(res, eigs)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)
            mat = utils.qid_from_idx(idx_u, np.prod(res))*mat

        elif field_row == ("velocityy",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)
            mat = utils.qid_from_idx(idx_v, np.prod(res))*mat

        elif field_row == ("velocityz",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)
            mat = utils.qid_from_idx(idx_w, np.prod(res))*mat

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)

        elif field_row == ("pressure",""):
            mat = c3d.zblk(res[0], res[1], res[2], 1, 1, 1, bc)

        return mat

    def linear_block(self, res, eq_params, eigs, bcs, field_row, field_col, restriction = None):
        """Create matrix block linear operator"""

        Ra = eq_params['rayleigh']
        Ta = eq_params['taylor']
        T = Ta**0.5

        xscale = eq_params['scale1d']
        yscale = eq_params['scale2d']
        zscale = eq_params['scale3d']

        idx_u, idx_v, idx_w, idx_p = self.zero_blocks(res, eigs, restriction = restriction)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_col)
        if field_row == ("velocityx",""):
            if field_col == ("velocityx",""):
                mat = c3d.i2j2k2lapl(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_u, np.prod(res))*mat*utils.qid_from_idx(idx_u, np.prod(res))
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + utils.id_from_idx(idx_u, np.prod(res))

            elif field_col == ("velocityy",""):
                mat = c3d.i2j2k2(res[0], res[1], res[2], bc, T, restriction = restriction)
                mat = utils.qid_from_idx(idx_u, np.prod(res))*mat*utils.qid_from_idx(idx_v, np.prod(res))

            elif field_col == ("velocityz",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = c3d.i2j2k2d1(res[0], res[1], res[2], bc, -1.0, xscale = xscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_u, np.prod(res))*mat*utils.qid_from_idx(idx_p, np.prod(res))

        elif field_row == ("velocityy",""):
            if field_col == ("velocityx",""):
                mat = c3d.i2j2k2(res[0], res[1], res[2], bc, -T, restriction = restriction)
                mat = utils.qid_from_idx(idx_v, np.prod(res))*mat*utils.qid_from_idx(idx_u, np.prod(res))

            elif field_col == ("velocityy",""):
                mat = c3d.i2j2k2lapl(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_v, np.prod(res))*mat*utils.qid_from_idx(idx_v, np.prod(res))
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + utils.id_from_idx(idx_v, np.prod(res))

            elif field_col == ("velocityz",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("pressure",""):
                mat = c3d.i2j2k2e1(res[0], res[1], res[2], bc, -1.0, yscale = yscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_v, np.prod(res))*mat*utils.qid_from_idx(idx_p, np.prod(res))

        elif field_row == ("velocityz",""):
            if field_col == ("velocityx",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("velocityz",""):
                mat = c3d.i2j2k2lapl(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_w, np.prod(res))*mat*utils.qid_from_idx(idx_w, np.prod(res))
                if bcs["bcType"] == self.SOLVER_HAS_BC:
                    mat = mat + utils.id_from_idx(idx_w, np.prod(res))

            elif field_col == ("temperature",""):
                mat = c3d.i2j2k2(res[0], res[1], res[2], bc, Ra, restriction = restriction)
                mat = utils.qid_from_idx(idx_w, np.prod(res))*mat

            elif field_col == ("pressure",""):
                mat = c3d.i2j2k2f1(res[0], res[1], res[2], bc, -1.0, zscale = zscale, restriction = restriction)
                mat = utils.qid_from_idx(idx_w, np.prod(res))*mat*utils.qid_from_idx(idx_p, np.prod(res))

        elif field_row == ("temperature",""):
            if field_col == ("velocityx",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("velocityy",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("velocityz",""):
                if self.linearize:
                    mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)
                    mat = mat*utils.qid_from_idx(idx_w, np.prod(res))
                else:
                    mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

            elif field_col == ("temperature",""):
                mat = c3d.i2j2k2lapl(res[0], res[1], res[2], bc, xscale = xscale, yscale = yscale, zscale = zscale, restriction = restriction)

            elif field_col == ("pressure",""):
                mat = c3d.zblk(res[0], res[1], res[2], 2, 2, 2, bc)

        elif field_row == ("pressure",""):
            if bcs["bcType"] == self.SOLVER_NO_TAU:
                mat = c3d.zblk(res[0], res[1], res[2], 1, 1, 1, no_bc())
            else:
                if field_col == ("velocityx",""):
                    bc['x']['cr'] = 1
                    bc['x']['rt'] = 1
                    bc['x']['zb'] = 1
                    bc['y']['cr'] = 1
                    bc['y']['rt'] = 1
                    bc['y']['zb'] = 1
                    bc['z']['cr'] = 1
                    bc['z']['rt'] = 1
                    bc['z']['zb'] = 1
                    mat = c3d.i1j1k1d1(res[0]+1, res[1]+1, res[2]+1, bc, xscale = xscale, restriction = restriction)
                    mat = utils.qid_from_idx(idx_p, np.prod(res))*mat*utils.qid_from_idx(idx_u, np.prod(res))

                elif field_col == ("velocityy",""):
                    bc['x']['cr'] = 1
                    bc['x']['rt'] = 1
                    bc['x']['zb'] = 1
                    bc['y']['cr'] = 1
                    bc['y']['rt'] = 1
                    bc['y']['zb'] = 1
                    bc['z']['cr'] = 1
                    bc['z']['rt'] = 1
                    bc['z']['zb'] = 1
                    mat = c3d.i1j1k1e1(res[0]+1, res[1]+1, res[2]+1, bc, yscale = yscale, restriction = restriction)
                    mat = utils.qid_from_idx(idx_p, np.prod(res))*mat*utils.qid_from_idx(idx_v, np.prod(res))

                elif field_col == ("velocityz",""):
                    bc['x']['cr'] = 1
                    bc['x']['rt'] = 1
                    bc['x']['zb'] = 1
                    bc['y']['cr'] = 1
                    bc['y']['rt'] = 1
                    bc['y']['zb'] = 1
                    bc['z']['cr'] = 1
                    bc['z']['rt'] = 1
                    bc['z']['zb'] = 1
                    mat = c3d.i1j1k1f1(res[0]+1, res[1]+1, res[2]+1, bc, zscale = zscale, restriction = restriction)
                    mat = utils.qid_from_idx(idx_p, np.prod(res))*mat*utils.qid_from_idx(idx_w, np.prod(res))

                elif field_col == ("temperature",""):
                    mat = c3d.zblk(res[0], res[1], res[2], 1, 1, 1, bc)

                elif field_col == ("pressure",""):
                    mat = c3d.zblk(res[0], res[1], res[2], 1, 1, 1, bc)
                    mat = mat + utils.id_from_idx(idx_p, np.prod(res))

        return mat

    def time_block(self, res, eq_params, eigs, bcs, field_row, restriction = None):
        """Create matrix block of time operator"""

        Pr = eq_params['prandtl']

        idx_u, idx_v, idx_w, idx_p = self.zero_blocks(res, eigs, restriction = restriction)

        bc = self.convert_bc(eq_params,eigs,bcs,field_row,field_row)
        if field_row == ("velocityx",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, 1.0/Pr, restriction = restriction)
            S = utils.qid_from_idx(idx_u, np.prod(res))
            mat = S*mat*S

        elif field_row == ("velocityy",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, 1.0/Pr, restriction = restriction)
            S = utils.qid_from_idx(idx_v, np.prod(res))
            mat = S*mat*S

        elif field_row == ("velocityz",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, 1.0/Pr, restriction = restriction)
            S = utils.qid_from_idx(idx_w, np.prod(res))
            mat = S*mat*S

        elif field_row == ("temperature",""):
            mat = c3d.i2j2k2(res[0], res[1], res[2], bc, restriction = restriction)

        elif field_row == ("pressure",""):
            mat = c3d.zblk(res[0], res[1], res[2], 1, 1, 1, bc)

        return mat

    def zero_blocks(self, res, eigs, restriction = None):
        """Build restriction matrices"""
    
        # U:
        idx_u = utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], res[1]-1), utils.qidx(res[2], 0), utils.qidx(res[0], 0))
        idx_u = np.union1d(idx_u, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], res[2]-1), utils.qidx(res[0], 0)))

        # V:
        idx_v = utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], 0), utils.qidx(res[0], res[0]-1))
        idx_v = np.union1d(idx_v, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], res[2]-1), utils.qidx(res[0], 0)))

        # W:
        idx_w = utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], 0), utils.qidx(res[0], res[0]-1))
        idx_w = np.union1d(idx_w,  utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], res[1]-1), utils.qidx(res[2], 0), utils.qidx(res[0], 0)))

        # Pressure: T_iNN, T_NjN, T_NNk
        idx_p = utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], 0), utils.qidx(res[0], res[0]-1))
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], res[1]-1), utils.qidx(res[2], 0), utils.qidx(res[0], 0)))
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], res[2]-1), utils.qidx(res[0], 0)))
        # Pressure: T_{N-2:N,N-2:N,N-2:N}
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], 0), utils.qidx(res[2], res[2]-3), utils.qidx(res[0], res[0]-3)))
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], res[1]-3), utils.qidx(res[2], res[2]-3), utils.qidx(res[0], 0)))
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.qidx(res[1], res[1]-3), utils.qidx(res[2], 0), utils.qidx(res[0], res[0]-3)))
        # Pressure: T_000
        idx_p = np.union1d(idx_p, utils.idx_kron_3d(res[1], res[2], res[0], utils.sidx(res[1], res[1]-1), utils.sidx(res[2], res[2]-1), utils.sidx(res[0], res[0]-1)))

        return (idx_u, idx_v, idx_w, idx_p)
