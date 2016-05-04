"""Module provides the functionalities to view linear stability solutions"""

from __future__ import division
from __future__ import unicode_literals

import numpy as np

abbr = {'prandtl':'Pr', 'rayleigh':'Ra', 'taylor':'Ta'}

def viewOperators(A, B = None, show = True, save = False):
    """Spy and/or write the operators to MatrixMarket file"""

    if save:
        import scipy.io as sciio
        sciio.mmwrite("matrix_A.mtx", A)
        sciio.mmwrite("matrix_B.mtx", B)

    # Spy the A (and B) operators
    if show:
        import matplotlib.pylab as pl

        if B is not None:
            pl.subplot(1,2,1)
        pl.spy(A, markersize=5, marker = '.', markeredgecolor = 'b')
        pl.tick_params(axis='x', labelsize=30)
        pl.tick_params(axis='y', labelsize=30)

        if B is not None:
            pl.subplot(1,2,2)
            pl.spy(B, markersize=5, marker = '.', markeredgecolor = 'b')
            pl.tick_params(axis='x', labelsize=30)
            pl.tick_params(axis='y', labelsize=30)

        pl.show()
        pl.clf()

def viewSpectra(fields, show = True, save = False, fid = None, max_cols = 3, subplot = False):
    """Plot spectra of eigenvectors"""

    if show or save:
        import matplotlib.pylab as pl
        # Plot spectra
        rows = np.ceil(len(fields)/max_cols)
        cols = min(max_cols, len(fields))
        for i,df in enumerate(fields.items()):
            if subplot:
                pl.subplot(rows,cols,i+1)
            pl.semilogy(np.abs(df[1]), 'k')
            #pl.semilogy(np.abs(df[1].real), ':')
            #pl.semilogy(np.abs(df[1].imag), ':')
            title = df[0][0]
            if df[0][1] != "":
                title = title + ', ' + df[0][1]
            pl.title(title)
            if save and not subplot:
                pl.tight_layout()
                fname = "spectra_" + title.replace(', ', '_')
                if fid is not None:
                    fname = fname + "_" + fid
                fname = fname + ".pdf"
                pl.savefig(fname, bbox_inches='tight', dpi=200)
            if show:
                pl.show()
            pl.clf()
        if subplot:
            pl.tight_layout()

            if save:
                fname = "spectra"
                if fid is not None:
                    fname = fname + "_" + fid
                fname = fname + ".pdf"
                pl.savefig(fname, bbox_inches='tight', dpi=200)

            if show:
                pl.show()
        pl.clf()

def viewPhysical(fields, geometry, res, eigs, eq_params, show = True, save = False, fid = None, max_cols = 3):
    """Plot eigenvectors in physical space"""

    if geometry == 'c1d':
        import geomhdiscc.transform.cartesian as transf
        nD = 1

        if ("pressure","") in fields:
            addContinuityC1D(fields, res, eigs, eq_params)
            viewSpectra(fields, show = show, save = save, fid=fid, max_cols = max_cols)

    elif geometry == 's1d':
        import geomhdiscc.transform.shell as transf
        nD = 1

    elif geometry == 'b1d':
        import geomhdiscc.transform.sphere as transf
        nD = 1

    elif geometry == 'c2d':
        import geomhdiscc.transform.cartesian as transf
        nD = 2

        if ("pressure","") in fields:
            addContinuityC2D(fields, res, eigs, eq_params)
            viewSpectra(fields, show = show, save = save, fid=fid, max_cols = max_cols)

    elif geometry == 'shell':
        import geomhdiscc.transform.shell as transf
        nD = 2

    elif geometry == "sphere":
        import geomhdiscc.transform.sphere as transf
        nD = 2

    elif geometry == "annulus":
        import geomhdiscc.transform.annulus as transf
        nD = 2

    elif geometry == "cylinder":
        import geomhdiscc.transform.cylinder_worland as transf
        nD = 2

    elif geometry == 'c3d':
        import geomhdiscc.transform.cartesian as transf
        nD = 3

        if ("pressure","") in fields:
            addContinuityC3D(fields, res, eq_params)
            viewSpectra(fields, show = show, save = save, fid=fid, max_cols = max_cols)

    if nD == 1:
        viewPhysical1D(fields, geometry, res, eigs, eq_params, transf, show = show, save = save, fid = fid, max_cols = max_cols)
    elif nD == 2:
        viewPhysical2D(fields, geometry, res, eigs, eq_params, transf, show = show, save = save, fid = fid, max_cols = max_cols)
    elif nD == 3:
        viewPhysical3D(fields, geometry, res, eigs, eq_params, transf, show = show, save = save, fid = fid, max_cols = max_cols)

def viewPhysical1D(specs, geometry, res, eigs, eq_params, transf, show = True, save = False, fid = None, max_cols = 3):
    """View 1D physical data"""

    sol_profile = dict()
    if geometry == 'c1d':
        viz_res = (res[0], )
        prof_opt = ()

    elif geometry == 's1d':
        import geomhdiscc.geometry.spherical.shell_radius as geo
        a, b = geo.linear_r2x(eq_params['ro'], eq_params['rratio'])
        viz_res = (res[0], a, b)
        prof_opt = ()

    elif geometry == 'b1d':
        viz_res = (res[0],)
        prof_opt = (int(eigs[0])%2,)

    for k,f in specs.items():
        sol_profile[k] = transf.toprofile(f, *prof_opt)

    if show or save:
        # Plot physical field
        grid = transf.grid_1d(*viz_res)
        viewProfile(sol_profile, grid, show = show, save = save, fid = fid, max_cols = max_cols)

def viewPhysical2D(specs, geometry, res, eigs, eq_params, transf, show = True, save = False, save_fast_profile = True, save_slow_profile = True, fid = None, max_cols = 3, slice_ratio = 4):
    """View 2D physical data"""

    sol_slice = dict()
    if geometry == 'c2d':
        res_1d = (res[0], )
        res_2d = (res[-1], )

    elif geometry == 'shell':
        import geomhdiscc.geometry.spherical.shell_radius as geo
        a, b = geo.linear_r2x(eq_params['ro'], eq_params['rratio'])
        res_1d = (res[0], a, b)
        res_2d = (res[1]-1, int(eigs[0]))

    elif geometry == 'sphere':
        res_1d = (res[0],)
        res_2d = (res[1]-1, int(eigs[0]))

    elif geometry == 'annulus':
        import geomhdiscc.geometry.cylindrical.annulus_radius as geo
        a, b = geo.linear_r2x(eq_params['ro'], eq_params['rratio'])
        res_1d = (res[0], a, b)
        res_2d = (res[-1],)

    elif geometry == 'cylinder':
        res_1d = (res[0],)
        res_2d = (int(eigs[0]), res[-1])

    viz_res = res_1d + res_2d

    for k,f in specs.items():
        sol_slice[k] = transf.toslice(f, res_1d[0], *res_2d) 

    if show or save:
        # Plot physical field as solution slice
        grid = transf.grid_2d(*viz_res)
        sfid = "solution"
        if fid is not None:
            sfid = sfid + "_" + fid
        viewSlice(sol_slice, grid, show = show, save = save, fid = sfid, max_cols = max_cols)

    if show or save:
        # Plot physical field on profile along fast direction
        grid_fast = transf.grid_fast(*res_1d)
        prof_fast = dict()
        for k,f in sol_slice.items():
            prof_fast[k] = f[np.floor(f.shape[0]/slice_ratio),:]

        pfid = "fast"
        if fid is not None:
            pfid = pfid + "_" + fid
        viewProfile(prof_fast, grid_fast, show = show, save = save, fid = pfid, max_cols = max_cols)

        if save_fast_profile:
            saveProfileData(prof_fast, grid_fast, fid = pfid)

        # Plot profile extruded along periodic direction
        grid_per = transf.grid_fast_per(*res_1d, m = np.ceil(eigs[0]))
        phi = eigs[0]*transf.eqgrid(np.ceil(eigs[0]))
        viewPeriodic(prof_fast, grid_per, phi, show = show, save = save, fid = pfid, max_cols = max_cols)

    if show or save:
        # Plot physical field on profile along slow direction
        grid_slow = transf.grid_slow(*res_2d)
        prof_slow = dict()
        for k,f in sol_slice.items():
            prof_slow[k] = f[:,np.floor(f.shape[1]/slice_ratio)]

        pfid = "slow"
        if fid is not None:
            pfid = pfid + "_" + fid
        viewProfile(prof_slow, grid_slow, show = show, save = save, fid = pfid, max_cols = max_cols)

        if save_slow_profile:
            saveProfileData(prof_slow, grid_slow, fid = pfid)

        # Plot profile extruded along periodic direction
        grid_per = transf.grid_slow_per(*res_2d, m = np.ceil(eigs[0]))
        phi = eigs[0]*transf.eqgrid(np.ceil(eigs[0]))
        viewPeriodic(prof_slow, grid_per, phi, show = show, save = save, fid = pfid, max_cols = max_cols)

def viewPhysical3D(specs, geometry, res, eigs, eq_params, transf, show = True, save = False, fid = None, max_cols = 3, slice_ratio = 2):
    """View 3D physical data"""

    sol_volume = dict()
    if geometry == 'c3d':
        res_1d = (res[0], )
        res_2d = (res[1], )
        res_3d = (res[2], )

    viz_res = res_1d + res_2d + res_3d

    for k,f in specs.items():
        sol_volume[k] = transf.tovolume(f, *viz_res) 

    if show or save:
        sol_slice = dict()
        # Plot physical field as slow slice
        slice_res = res_1d + res_2d
        grid = transf.grid_2d(*slice_res)
        for k,f in sol_volume.items():
            sol_slice[k] = f[:,:,np.floor(f.shape[2]/slice_ratio)].T
        sfid = "slow"
        if fid is not None:
            sfid = sfid + "_" + fid
        viewSlice(sol_slice, grid, show = show, save = save, fid = sfid, max_cols = max_cols)

        sol_slice = dict()
        # Plot physical field as medium slice
        slice_res = res_1d + res_3d
        grid = transf.grid_2d(*slice_res)
        for k,f in sol_volume.items():
            sol_slice[k] = f[:,np.floor(f.shape[1]/slice_ratio),:].T 
        sfid = "medium"
        if fid is not None:
            sfid = sfid + "_" + fid
        viewSlice(sol_slice, grid, show = show, save = save, fid = sfid, max_cols = max_cols)

        sol_slice = dict()
        # Plot physical field as fast slice
        slice_res = res_2d + res_3d
        grid = transf.grid_2d(*slice_res)
        for k,f in sol_volume.items():
            sol_slice[k] = f[np.floor(f.shape[0]/slice_ratio),:,:].T
        sfid = "fast"
        if fid is not None:
            sfid = sfid + "_" + fid
        viewSlice(sol_slice, grid, show = show, save = save, fid = sfid, max_cols = max_cols)

def saveProfileData(fields, grid, fid = None):
    """Save profile data to ASCII file"""

    fbase = "profile"
    if fid is not None:
        fbase = fbase + "_" + fid

    for k,f in fields.items():
        fname = fbase + "_" + k[0]
        if k[1] != "":
            fname = fname + "_" + k[1]
        np.savetxt(fname + ".dat", np.array([grid, f.real, f.imag]).transpose())

def viewProfile(fields, grid, fid = None, show = True, save = False, max_cols = 3):
    """View a profile"""

    import matplotlib.pylab as pl
    rows = np.ceil(len(fields)/max_cols)
    cols = min(max_cols, len(fields))
    for i,df in enumerate(fields.items()):
        pl.subplot(rows,cols,i+1)
        pl.plot(grid, df[1].real, 'b-')
        pl.plot(grid, df[1].imag, 'r-')
        title = df[0][0]
        if df[0][1] != "":
            title = title + ', ' + df[0][1]
        pl.title(title)
    pl.tight_layout()
    if save:
        fname = "profile"
        if fid is not None:
            fname = fname + "_" + fid
        fname = fname + ".pdf"
        pl.savefig(fname, bbox_inches='tight', dpi=200)
    if show:
        pl.show()
    pl.clf()

def viewSlice(fields, grid, fid = None, show = True, save = False, max_cols = 3, subplot = False):
    """View a slice"""

    import matplotlib as mpl
    import matplotlib.pylab as pl
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # Plot physical field in meridional slice
    rows = np.ceil(len(fields)/max_cols)
    cols = min(max_cols, len(fields))
    mycm = cm.bwr
    for i,df in enumerate(fields.items()):
        vmax = np.max(np.abs(df[1].real))
        if subplot:
            pl.subplot(rows,cols,i+1, aspect = 'equal', axisbg = 'black')
        CS = pl.contourf(grid[0], grid[1], df[1].real, 30, cmap = mycm, vmax = vmax, vmin = -vmax)
        title = df[0][0]
        if df[0][1] != "":
            title = title + ', ' + df[0][1]
        pl.title(title)
        divider = make_axes_locatable(pl.gca())
        cax = divider.append_axes("right", "5%", pad="3%")
        pl.colorbar(CS, cax=cax)
        if not subplot:
            pl.tight_layout()
            fname = "slice_" + title.replace(', ', '_')
            if fid is not None:
                fname = fname + "_" + fid
            fname = fname + ".pdf"
            pl.savefig(fname, bbox_inches='tight', dpi=200)
            if show:
                pl.show()
            pl.clf()
    if subplot:
        pl.tight_layout()
        if save:
            fname = "slice"
            if fid is not None:
                fname = fname + "_" + fid
            fname = fname + ".pdf"
            pl.savefig(fname, bbox_inches='tight', dpi=200)
        if show:
            pl.show()
    pl.clf()

def viewPeriodic(fields, grid, grid_per, fid = None, show = True, save = False, max_cols = 3):
    """Extend field profile along periodci direction"""

    sol_per = dict()
    for k,f in fields.items():
        sol_per[k] = np.outer(np.cos(grid_per),f.real) - np.outer(np.sin(grid_per), f.imag)

    viewSlice(sol_per, grid, fid = fid, show = show, save = save, max_cols = max_cols)

def addContinuityC1D(fields, res, eigs, eq_params):
    """Add continuity fields for cartesian 1D"""

    import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
    from geomhdiscc.geometry.cartesian.cartesian_boundary_1d import no_bc
    if len(res) == 3:
        k1 = eigs[0]
        k2 = eigs[1]
    else:
        k1 = eigs[0]
        k2 = 0 
    mat = c1d.i2lapl(res[0], k1, k2, no_bc(), cscale = eq_params['scale1d'])
    f = fields[("pressure","")]
    fields[("pressure","elliptic")] = mat*f
    if ("temperature","") in fields:
        mat = c1d.i2d1(res[0], no_bc(), -eq_params['rayleigh'], cscale = eq_params['scale1d'])
        f = fields[("temperature","")]
        fields[("pressure","elliptic")] = fields[("pressure","elliptic")] + mat*f

    mat = c1d.i1d1(res[0], no_bc(), cscale = eq_params['scale1d'])
    cont = mat*fields[("velocity","z")]
    if ("velocity","x") in fields: 
        mat = c1d.i1(res[0], no_bc(), 1j*k1)
        cont = cont + mat*fields[("velocity","x")]
    if ("velocity","y") in fields: 
        mat = c1d.i1(res[0], no_bc(), 1j*k2)
        cont = cont + mat*fields[("velocity","y")]
    fields[("continuity","")] = cont

def addContinuityC2D(fields, res, eigs, eq_params):
    """Add continuity fields for cartesian 2D"""

    import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
    from geomhdiscc.geometry.cartesian.cartesian_boundary_2d import no_bc
    xscale = eq_params['scale1d']
    if len(res) == 3:
        zscale = eq_params['scale3d']
    else:
        zscale = eq_params['scale2d']

    mat = c2d.i2j2lapl(res[0], res[-1], eigs[0], no_bc(), xscale = xscale, zscale = zscale)
    f = fields[("pressure","")]
    fields[("pressure","elliptic")] = mat*f
    if ("temperature","") in fields:
        mat = c2d.i2j2e1(res[0], res[-1], no_bc(), -eq_params['rayleigh'], zscale = zscale)
        f = fields[("temperature","")]
        fields[("pressure","elliptic")] = fields[("pressure","elliptic")] + mat*f

    mat = c2d.i1j1d1(res[0], res[-1], no_bc(), xscale = xscale)
    cont = mat*fields[("velocity","x")]
    mat = c2d.i1j1e1(res[0], res[-1], no_bc(), zscale = zscale)
    cont = cont +  mat*fields[("velocity","z")]
    if ("velocity","y") in fields: 
        mat = c2d.i1j1(res[0], res[-1], no_bc(), 1j*eigs[0])
        cont = cont + mat*fields[("velocity","y")]
    fields[("continuity","")] = cont

def addContinuityC3D(fields, res, eq_params):
    """Add continuity fields for cartesian 3D"""

    import geomhdiscc.geometry.cartesian.cartesian_3d as c3d
    from geomhdiscc.geometry.cartesian.cartesian_boundary_3d import no_bc
    mat = c3d.i2j2k2lapl(res[0], res[1], res[2], no_bc(), xscale = eq_params['scale1d'], yscale = eq_params['scale2d'], zscale = eq_params['scale3d'])
    f = fields[("pressure","")]
    fields[("pressure","elliptic")] = mat*f
    if ("temperature","") in fields:
        mat = c3d.i2j2k2f1(res[0], res[1], res[2], no_bc(), -eq_params['rayleigh'], zscale = eq_params['scale3d'])
        f = fields[("temperature","")]
        fields[("pressure","elliptic")] = fields[("pressure","elliptic")] + mat*f

    mat = c3d.i1j1k1d1(res[0], res[1], res[2], no_bc(), xscale = eq_params['scale1d'])
    cont = mat*fields[("velocity","x")]
    mat = c3d.i1j1k1e1(res[0], res[1], res[2], no_bc(), yscale = eq_params['scale3d'])
    cont = cont +  mat*fields[("velocity","y")]
    mat = c3d.i1j1k1f1(res[0], res[1], res[2], no_bc(), zscale = eq_params['scale3d'])
    cont = cont + mat*fields[("velocity","z")]
    fields[("continuity","")] = cont
