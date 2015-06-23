"""Module provides the functor required to trace a marginal curve"""

from __future__ import division
from __future__ import unicode_literals

import copy
import numpy as np
import scipy.optimize as optimize
import functools

import geomhdiscc.linear_stability.solver_slepc as solver_mod
import geomhdiscc.linear_stability.io as io
from geomhdiscc.linear_stability.solver_slepc import Print
from mpi4py import MPI


class MarginalCurve:
    """The marginal curve"""

    def __init__(self, gevp_opts, mode, rtol = 1e-7, evp_tol = 1e-8):
        """Initialize the marginal curve"""

        # Open file for IO and write header
        if MPI.COMM_WORLD.Get_rank() == 0:
            self.out = open('marginal_curve.dat', 'a', 0)
        else:
            self.out = None
        io.write_header(self.out, gevp_opts['model'].__class__.__name__, len(gevp_opts['res']), len(gevp_opts['eigs']), gevp_opts['eq_params'])

        self.point = MarginalPoint(gevp_opts, mode, rtol = rtol, out_file = self.out, evp_tol = evp_tol)
        self.rtol = rtol
        self.guess = None

    def trace(self, ks, guess = None, initial_guess = 1.0):
        """Trace the marginal curve"""
        
        Print("Tracing marginal curve")
        Print("----------------------")

        data_k = np.zeros(len(ks))
        data_Ra = np.zeros(len(ks))
        data_freq = np.zeros(len(ks))
        Rac = guess
        for i, k in enumerate(ks):
            if i > 3:
                Rac = p(k)
            Rac, freq = self.point(k, guess = Rac, initial_guess = initial_guess)
            data_k[i] = k
            data_Ra[i] = Rac
            data_freq[i] = freq
            if i > 2:
                lb = max(0, i-3)
                z = np.polyfit(data_k[lb:i+1], data_Ra[lb:i+1], min(i,3))
                p = np.poly1d(z)

        return (data_k, data_Ra, data_freq)

    def minimum(self, ks, Ras, xtol = 1e-4, only_int = False):
        """Compute the critical value (minimum of curve)"""

        # Open file for IO and write header
        if MPI.COMM_WORLD.Get_rank() == 0:
            min_out = open('marginal_minimum.dat', 'a', 0)
        else:
            min_out = None
        io.write_header(min_out, self.point._gevp.model.__class__.__name__, len(self.point._gevp.res), len(self.point._gevp.eigs), self.point._gevp.eq_params)
        
        Print("Finding minimum")
        Print("---------------")

        # Get bounds for minimum
        imin = np.argmin(Ras)
        a = ks[imin-1]
        b = ks[imin]
        c = ks[imin+1]
        self.guess = Ras[imin]
        if imin == 0 or imin == len(Ras)-1:
            Print("Minimum is not in provided range")
            return (None, None, None)
        else:
            if only_int:
                self.point(b, guess = Ras[imin])
                min_kc = b
                min_Rac = Ras[imin]
                min_fc = self.point.frequency()

            else:
                # Find minimum
                opt = optimize.minimize_scalar(self._tracker, [a, b, c], options = {'xtol':xtol})
                min_kc = opt.x
                min_Rac = opt.fun
                min_fc = self.point.frequency()

            Print("Minimumt at kc: {:g}, Rac: {:g}, fc: {:g}".format(min_kc, min_Rac, min_fc))
            io.write_results(min_out, self.point._gevp.res, self.point._gevp.wave(min_kc), self.point._gevp.eq_params, min_Rac, min_fc, self.point.mode, self.point.growthRate())
            return (min_kc, min_Rac, min_fc)

    def view(self, ks, Ras, fs, minimum = None, plot = True, save_pdf = False):
        """Plot the marginal curve"""

        if plot or save_pdf:
            import matplotlib.pylab as pl

            pl.subplot(1,2,1)
            pl.plot(ks, Ras, 'b')
            if minimum is not None:
                pl.plot(minimum[0], minimum[1], 'r+', markersize=14)
            pl.title('Critical Rayleigh number')
            pl.xlabel('Wavenumber')
            pl.ylabel('Ra')
            pl.tight_layout()

            pl.subplot(1,2,2)
            pl.plot(ks, fs, 'g')
            pl.title('Critical frequency')
            pl.xlabel('Wavenumber')
            pl.ylabel('frequency')
            if save_pdf:
                pl.savefig('marginal_curve.pdf', bbox_inches='tight', dpi = 200)
            if plot:
                pl.tight_layout()
                pl.show()
                pl.clf()

    def _tracker(self, k):

        return self.point(k, guess = self.guess)[0]

class MarginalPoint:
    """A point on the marginal curve"""

    def __init__(self, gevp_opts, mode, rtol = 1e-7, out_file = None, evp_tol = 1e-8):
        """Initialize the marginal point"""

        self.out = out_file
        gevp_opts['tol'] = evp_tol
        self._gevp = GEVP(**gevp_opts)
        self.mode = mode
        self.Rac = None
        self.fc = None
        self.a = None
        self.b = None
        self.rtol = rtol
        self.nev = max(mode+1, 1)
        self.best = None

    def __call__(self, k, guess = None, initial_guess = 1.0):
        """Compute marginal point in given interval [a, b]"""

        if guess is None:
            Ra = initial_guess
            factor = 2.0
        else:
            Ra = guess
            factor = 1.1

        Print("Computing marginal curve at k = {:g}".format(k))
        Print("\t- initial guess Ra = {:g}".format(Ra))

        # Set wavenumter
        self._gevp.setEigs(k)

        # Search interval (newton)
        self._gevp.best = None
        Print("\t- finding interval")
        for i in range(0,16):
            self.findInterval(Ra)
            if self.a != None and self.b != None:
                break
            elif self.best is not None and self.best[2]:
                Ra = self.best[1]
                self.best = (self.best[0], self.best[1], False)
            elif self.best is not None:
                if self.best[0] > 0:
                    Ra = self.best[1]/(1.1*factor)
                else:
                    Ra = self.best[1]*(1.1*factor)
                self.best = None
            else:
                Ra = factor*Ra

        if self.a is None or self.b is None:
            self.Rac = float('NaN')
            self.rc = float('NaN')
            self.a = None
            self.b = None
        else:
            # Refine zero with Brent's method
            Print("\t- Brent's method: [{:g}, {:g}]".format(self.a,self.b))
            Rac = optimize.brentq(self._tracker, self.a, self.b, rtol = self.rtol)
            self.a = None
            self.b = None
            self.Rac = Rac
            self.fc = self.frequency()

        # Write results to file and show on console output
        io.write_results(self.out, self._gevp.res, self._gevp.eigs, self._gevp.eq_params, self.Rac, self.fc, self.mode, self.growthRate())
        Print("\t- Ra = {:g}, freq = {:g}, conv = {:g}".format(self.Rac, self.fc, self.growthRate()))

        return (self.Rac, self.fc)

    def findInterval(self, guess):
        """Find the sign changing interval"""

        try:
            #Rac = optimize.newton(self._signTracker, guess)
            Rac = solver_mod.newton(self._signTracker, guess, step = 1e-4)
        except (solver_mod.NewtonDoneError, solver_mod.NewtonNoneError):
            pass
        except:
            raise

    def eigenvector(self):
        """Compute the eigenvectors"""

        self._gevp.solve(self.Rac, self.nev, with_vectors = True)

        return self._gevp.evp_vec[:,self.mode]

    def eigenvalue(self):
        """Compute the eigenvalues"""

        self._gevp.solve(self.Rac, self.nev)

        return self._gevp.evp_lmb[self.mode]

    def growthRate(self):
        """Get the growth rate"""

        self._gevp.solve(self.Rac, self.nev)

        return self._gevp.evp_lmb[self.mode].real

    def frequency(self):
        """Compute the frequency rate"""

        self._gevp.solve(self.Rac, self.nev, use_vector = self.mode)

        return self._gevp.evp_lmb[self.mode].imag

    def _tracker(self, Ra):
        """Track the growth rate to find critical value"""

        for i in range(0, 10):
            self._gevp.solve(Ra, self.nev, with_vectors = True, use_vector = self.mode)
            if self._gevp.evp_lmb is not None:
                break

        growth = self._gevp.evp_lmb[self.mode].real
        if self.best is None:
            self.best = (abs(growth), Ra, False)
        elif abs(growth) < self.best[0]:
            self.best = (abs(growth), Ra, True)

        return growth

    def _signTracker(self, Ra):
        """Track the growth rate to find critical value"""

        # Get change in rayleigh number
        if self._gevp.eq_params['rayleigh'] is None:
            dRa = Ra
        else:
            dRa = Ra - self._gevp.eq_params['rayleigh']

        # Solve
        self._gevp.solve(Ra, self.nev)

        if self._gevp.evp_lmb is None:
            raise solver_mod.NewtonNoneError

        growth = self._gevp.evp_lmb[self.mode].real
        if growth > 0 and growth < 100:
            self.b = Ra
        elif growth < 0 and growth > -100:
            self.a = Ra

        if self.a is not None and self.b is not None:
            raise solver_mod.NewtonDoneError

        # Check for change of sign (require small change in rayleigh and small growth rate)
        if abs(dRa/Ra) < 0.01 and abs(self._gevp.evp_lmb[self.mode].real) < 0.1:
            self._gevp.solve(Ra+dRa, self.nev)
            if self._gevp.evp_lmb is not None and growth*self._gevp.evp_lmb[self.mode].real < 0:
                self.a = min(Ra, Ra+dRa)
                self.b = max(Ra, Ra+dRa)
                raise solver_mod.NewtonDoneError

        if self.best is None:
            self.best = (abs(growth), Ra, False)
        elif abs(growth) < self.best[0]:
            self.best = (abs(growth), Ra, True)

        return growth

class GEVP:
    """Class to represent a marginal point on a marginal curve"""
    
    def __init__(self, model, res, eq_params, eigs, bcs, wave = None, tol = 1e-8, ellipse_radius = None, fixed_shift = False):
        """Initialize the marginal point variables"""

        self.model = copy.copy(model)
        self.res = res
        self.eq_params = eq_params.copy()
        self.eq_params['rayleigh'] = None
        self.eigs = eigs
        self.bcs = bcs.copy()
        self.no_bcs = bcs.copy()
        self.no_bcs['bcType'] = model.SOLVER_NO_TAU
        self.fields = self.model.stability_fields()
        self.evp_lmb = None
        self.evp_vec = None
        self.changed = True
        self.solver = solver_mod.GEVPSolver(tol = tol, ellipse_radius = ellipse_radius, fixed_shift = fixed_shift)
        if wave is None:
            self.wave = self.defaultWave
        else:
            self.wave = wave

    def setupProblem(self, Ra):
        """Build A and B matrices for the given rayleigh number"""

        self.eq_params['rayleigh'] = Ra

        opA = functools.partial(self.model.implicit_linear, self.res, self.eq_params, self.eigs, self.bcs, self.fields)
        opB = functools.partial(self.model.time, self.res, self.eq_params, self.eigs, self.no_bcs, self.fields)
        sizes = self.model.stability_sizes(self.res, self.eigs)

        return (opA, opB, sizes)
   
    def setEigs(self, k):
       """Set the wave number"""

       self.eigs = self.wave(k)
       self.changed = True

    def solve(self, Ra, nev, with_vectors = False, use_vector = None):
        """Solve the GEVP and store eigenvalues"""

        if self.changed or Ra != self.eq_params['rayleigh'] or (with_vectors and self.evp_vec is None):
            self.changed = False
            problem = self.setupProblem(Ra)

            if use_vector is not None and self.evp_vec is not None:
                initial_vector = self.evp_vec[:,use_vector]
            else:
                initial_vector = None

            if with_vectors:
                self.evp_lmb, self.evp_vec = self.solver.eigenpairs(problem, nev, initial_vector = initial_vector)
            else:
                self.evp_lmb = self.solver.eigenvalues(problem, nev, initial_vector = initial_vector)
                self.evp_vec = None
            Print("\t\t (Ra = {:.14g}, ev = ".format(Ra) + str(self.evp_lmb) + ")")

    def viewOperators(self, Ra, spy = True, write_mtx = True):
        """Spy and/or write the operators to MatrixMarket file"""

        # Create operators
        if spy or write_mtx:
            opA, opB, sizes = self.setupProblem(Ra)
            restrict = self.solver.restrict_operators(sizes)
            A = opA(restriction = restrict)
            B = opB(restriction = restrict)

        if write_mtx:
            import scipy.io as sciio
            sciio.mmwrite("matrix_A.mtx", A)
            sciio.mmwrite("matrix_B.mtx", B)

        # Spy the two operators
        if spy:
            import matplotlib.pylab as pl
            pl.subplot(1,2,1)
            pl.spy(A, markersize=5, marker = '.', markeredgecolor = 'b')
            pl.tick_params(axis='x', labelsize=30)
            pl.tick_params(axis='y', labelsize=30)
            pl.subplot(1,2,2)
            pl.spy(B, markersize=5, marker = '.', markeredgecolor = 'b')
            pl.tick_params(axis='x', labelsize=30)
            pl.tick_params(axis='y', labelsize=30)
            pl.show()
            pl.clf()

    def viewSpectra(self, viz_mode, plot = True, save_pdf = False):
        """Plot the spectra of the eigenvectors"""

        if self.evp_lmb is not None:

            # Extra different fields
            start = 0
            stop = 0
            sol_spec = dict()
            for i,f in enumerate(self.fields):
                stop = stop + self.model.stability_sizes(self.res, self.eigs)[0][i]
                sol_spec[f] = self.evp_vec[start:stop, viz_mode]
                start = stop
            
            if plot or save_pdf:
                import matplotlib.pylab as pl
                Print("\nVisualizing spectra of mode: " + str(self.evp_lmb[viz_mode]))
                # Plot spectra
                rows = np.ceil(len(self.fields)/3)
                cols = min(3, len(self.fields))
                for i,f in enumerate(self.fields):
                    pl.subplot(rows,cols,i+1)
                    pl.semilogy(np.abs(sol_spec[f]))
                    title = f[0]
                    if f[1] != "":
                        title = title + ', ' + f[1]
                    pl.title(title)
                    pl.tight_layout()
                if save_pdf:
                    fname = "spectra_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                    pl.savefig(fname, bbox_inches='tight', dpi=200)
                if plot:
                    pl.show()
                pl.clf()

            return sol_spec

        else:
            return None

    def viewPhysical(self, viz_mode, geometry, plot = True, save_pdf = False, save_profile = True):
        """Plot the spectra of the eigenvectors"""

        if self.evp_lmb is not None:

            # Extra different fields
            sol_spec = self.viewSpectra(viz_mode, plot = False, save_pdf = False)
            if self.model.use_galerkin:
                sol_spec = self.apply_stencil(sol_spec)

            if geometry == 'c1d':
                import geomhdiscc.transform.cartesian as transf
                nD = 1

                if ("pressure","") in self.fields:
                    import geomhdiscc.geometry.cartesian.cartesian_1d as c1d
                    mat = c1d.i2lapl(self.res[0], self.eigs[0], self.eigs[1], {0:0}, cscale = 2.0)
                    field = sol_spec[("pressure","")]
                    sol_spec[("pressure","laplacian")] = mat*field
                    mat = c1d.i2d1(self.res[0], {0:0}, cscale = 2.0)
                    field = sol_spec[("temperature","")]
                    sol_spec[("pressure","laplacian")] = sol_spec[("pressure","laplacian")] - self.eq_params['rayleigh']*mat*field

                    mat = c1d.i1d1(self.res[0], {0:0}, cscale = 2.0)
                    cont = mat*sol_spec[("velocity","z")]
                    if ("velocity","x") in self.fields: 
                        mat = c1d.i1(self.res[0], {0:0}, 1j*self.eigs[0])
                        cont = cont + mat*sol_spec[("velocity","x")]
                    if ("velocity","y") in self.fields: 
                        mat = c1d.i1(self.res[0], {0:0}, 1j*self.eigs[1])
                        cont = cont + mat*sol_spec[("velocity","y")]
                    sol_spec[("continuity","")] = cont

            elif geometry == 's1d':
                import geomhdiscc.transform.shell as transf
                nD = 1
            elif geometry == 'b1d':
                import geomhdiscc.transform.sphere as transf
                nD = 1
            elif geometry == 'c2d':
                import geomhdiscc.transform.cartesian as transf
                nD = 2

                if ("pressure","") in self.fields:
                    import geomhdiscc.geometry.cartesian.cartesian_2d as c2d
                    mat = c2d.i2j2lapl(self.res[0], self.res[2], self.eigs[0], {'x':{0:0}, 'z':{0:0}}, xscale = self.eq_params['scale1d'], zscale = self.eq_params['scale3d'])
                    field = sol_spec[("pressure","")]
                    sol_spec[("pressure","laplacian")] = mat*field
                    mat = c2d.i2j2e1(self.res[0], self.res[2], {'x':{0:0}, 'z':{0:0}}, zscale = self.eq_params['scale3d'])
                    field = sol_spec[("temperature","")]
                    sol_spec[("pressure","laplacian")] = sol_spec[("pressure","laplacian")] - self.eq_params['rayleigh']*mat*field

                    mat = c2d.i1j1d1(self.res[0], self.res[2], {'x':{0:0}, 'z':{0:0}}, xscale = self.eq_params['scale1d'])
                    cont = mat*sol_spec[("velocity","x")]
                    mat = c2d.i1j1e1(self.res[0], self.res[2], {'x':{0:0}, 'z':{0:0}}, zscale = self.eq_params['scale3d'])
                    cont = cont +  mat*sol_spec[("velocity","z")]
                    if ("velocity","y") in self.fields: 
                        mat = c2d.i1j1(self.res[0], self.res[2], {'x':{0:0}, 'z':{0:0}}, 1j*self.eigs[0])
                        cont = cont + mat*sol_spec[("velocity","y")]
                    sol_spec[("continuity","")] = cont

            elif geometry == 'shell':
                import geomhdiscc.transform.shell as transf
                nD = 2
            elif geometry == "sphere":
                import geomhdiscc.transform.sphere as transf
                nD = 2
            elif geometry == "annulus":
                import geomhdiscc.transform.annulus as transf
                nD = 2

            # 1D data: plot physical profile
            sol_phys = dict()
            if nD == 1:
                if geometry == 'c1d':
                    viz_res = (self.res[0], )
                    prof_opt = ()

                elif geometry == 's1d':
                    import geomhdiscc.geometry.spherical.shell_radius as geo
                    a, b = geo.linear_r2x(self.eq_params['ro'], self.eq_params['rratio'])
                    viz_res = (self.res[0], a, b)
                    prof_opt = ()

                elif geometry == 'b1d':
                    viz_res = (self.res[0],)
                    prof_opt = (int(self.eigs[0])%2,)

                for k,f in sol_spec.items():
                    sol_phys[k] = transf.toprofile(f, *prof_opt)

                if plot or save_pdf:
                    import matplotlib.pylab as pl
                    Print("\nVisualizing physical data of mode: " + str(self.evp_lmb[viz_mode]))
                    # Plot physical field
                    grid = transf.grid_1d(*viz_res)
                    rows = np.ceil(len(sol_phys)/3)
                    cols = min(3, len(sol_phys))
                    for i,df in enumerate(sol_phys.items()):
                        pl.subplot(rows,cols,i+1)
                        pl.plot(grid, df[1].real)
                        pl.plot(grid, df[1].imag)
                        title = df[0][0]
                        if df[0][1] != "":
                            title = title + ', ' + df[0][1]
                        pl.title(title)
                        pl.tight_layout()
                    if save_pdf:
                        fname = "profile_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()

            # 2D data: plot physical contours on slice
            elif nD == 2:
                if geometry == 'c2d':
                    res_1d = (self.res[0], )
                    res_2d = (self.res[2], )

                elif geometry == 'shell':
                    import geomhdiscc.geometry.spherical.shell_radius as geo
                    a, b = geo.linear_r2x(self.eq_params['ro'], self.eq_params['rratio'])
                    res_1d = (self.res[0], a, b)
                    res_2d = (self.res[1]-1, int(self.eigs[0]))

                elif geometry == 'sphere':
                    res_1d = (self.res[0],)
                    res_2d = (self.res[1]-1, int(self.eigs[0]))

                elif geometry == 'annulus':
                    import geomhdiscc.geometry.cylindrical.annulus_radius as geo
                    a, b = geo.linear_r2x(self.eq_params['ro'], self.eq_params['rratio'])
                    res_1d = (self.res[0], a, b)
                    res_2d = (self.res[2],)

                viz_res = res_1d + res_2d

                for k,f in sol_spec.items():
                    sol_phys[k] = transf.toslice(f, res_1d[0], *res_2d) 

                if plot or save_pdf:
                    import matplotlib as mpl
                    import matplotlib.pylab as pl
                    import matplotlib.cm as cm
                    from mpl_toolkits.axes_grid1 import make_axes_locatable
                    Print("\nVisualizing physical data of mode: " + str(self.evp_lmb[viz_mode]))
                    # Plot physical field in meridional slice
                    grid = transf.grid_2d(*viz_res)
                    rows = np.ceil(len(sol_phys)/3)
                    cols = min(3, len(sol_phys))
                    mycm = cm.bwr
                    for i,df in enumerate(sol_phys.items()):
                        vmax = np.max(np.abs(df[1].real))
                        pl.subplot(rows,cols,i+1, aspect = 'equal', axisbg = 'black')
                        CS = pl.contourf(grid[0], grid[1], df[1].real, 30, cmap = mycm, vmax = vmax, vmin = -vmax)
                        title = df[0][0]
                        if df[0][1] != "":
                            title = title + ', ' + df[0][1]
                        pl.title(title)
                        divider = make_axes_locatable(pl.gca())
                        cax = divider.append_axes("right", "5%", pad="3%")
                        pl.colorbar(CS, cax=cax)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "slice_meridional_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()

                if (geometry == 'shell' or geometry == 'sphere') and (plot or save_pdf or save_profile):
                    # Plot physical field in equatorial radial profile
                    grid_r = transf.rgrid(*res_1d)
                    sol_rrad = dict()
                    sol_irad = dict()
                    for i,f in enumerate(self.fields):
                        sol_rrad[f] = sol_phys[f].real[sol_phys[f].shape[0]//2,:]
                        sol_irad[f] = sol_phys[f].imag[sol_phys[f].shape[0]//2,:]
                        if f == ("temperature",""):
                            np.savetxt("profile_equatorial_Ta{:g}".format(self.eq_params['taylor']) + ".dat", np.array([grid_r, sol_rrad[f], sol_irad[f]]).transpose())

                if (geometry == 'shell' or geometry == 'sphere') and (plot or save_pdf):
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        pl.subplot(rows,cols,i+1)
                        pl.plot(grid_r, sol_rrad[f], 'b-')
                        pl.plot(grid_r, sol_irad[f], 'r-')
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "profile_equatorial_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()
                    # Plot physical field in equatorial slice
                    grid_eq = transf.grid_eq(*res_1d, m = int(self.eigs[0]))
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        phi = self.eigs[0]*transf.eqgrid(int(self.eigs[0]))
                        sol_eq = np.outer(np.cos(phi),sol_rrad[f]) - np.outer(np.sin(phi),sol_irad[f])
                        pl.subplot(rows,cols,i+1, aspect = 'equal', axisbg = 'black')
                        vmax = np.max(np.abs(sol_eq))
                        CS = pl.contourf(grid_eq[0], grid_eq[1], sol_eq, 30, cmap = mycm)
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                        divider = make_axes_locatable(pl.gca())
                        cax = divider.append_axes("right", "5%", pad="3%")
                        pl.colorbar(CS, cax=cax)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "slice_equatorial_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()

                if (geometry == 'c2d') and (plot or save_pdf or save_profile):
                    # Plot physical field in equatorial radial profile
                    grid_z = transf.grid(*res_2d)
                    sol_prof = dict()
                    for i,f in enumerate(self.fields):
                        sol_prof[f] = sol_phys[f][:,sol_phys[f].shape[1]//2]

                if (geometry == 'c2d') and (plot or save_pdf):
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        pl.subplot(rows,cols,i+1)
                        pl.plot(grid_z, sol_prof[f].real, 'b-')
                        pl.plot(grid_z, sol_prof[f].imag, 'r-')
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "profile_vertical_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()
                    # Plot physical field in vertical slice
                    grid_eq = transf.grid_eq(*res_2d, m = int(self.eigs[0]))
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        phi = self.eigs[0]*transf.eqgrid(int(self.eigs[0]))
                        sol_eq = np.outer(sol_prof[f].real,np.cos(phi)) - np.outer(sol_prof[f].imag, np.sin(phi))
                        pl.subplot(rows,cols,i+1, aspect = 'equal', axisbg = 'black')
                        vmax = np.max(np.abs(sol_eq))
                        CS = pl.contourf(grid_eq[0], grid_eq[1], sol_eq, 30, cmap = mycm)
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                        divider = make_axes_locatable(pl.gca())
                        cax = divider.append_axes("right", "5%", pad="3%")
                        pl.colorbar(CS, cax=cax)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "slice_equatorial_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()

                if (geometry == 'c2d') and (plot or save_pdf or save_profile):
                    # Plot physical field in equatorial radial profile
                    grid_z = transf.grid(*res_1d)
                    sol_prof = dict()
                    for i,f in enumerate(self.fields):
                        sol_prof[f] = sol_phys[f][sol_phys[f].shape[0]//2,:]

                if (geometry == 'c2d') and (plot or save_pdf):
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        pl.subplot(rows,cols,i+1)
                        pl.plot(grid_z, sol_prof[f].real, 'b-')
                        pl.plot(grid_z, sol_prof[f].imag, 'r-')
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "profile_vertical_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()
                    # Plot physical field in vertical slice
                    grid_eq = transf.grid_eq(*res_2d, m = int(self.eigs[0]))
                    rows = np.ceil(len(self.fields)/3)
                    cols = min(3, len(self.fields))
                    for i,f in enumerate(self.fields):
                        phi = self.eigs[0]*transf.eqgrid(int(self.eigs[0]))
                        sol_eq = np.outer(sol_prof[f].real,np.cos(phi)) - np.outer(sol_prof[f].imag, np.sin(phi))
                        pl.subplot(rows,cols,i+1, aspect = 'equal', axisbg = 'black')
                        vmax = np.max(np.abs(sol_eq))
                        CS = pl.contourf(grid_eq[0], grid_eq[1], sol_eq, 30, cmap = mycm)
                        title = f[0]
                        if f[1] != "":
                            title = title + ', ' + f[1]
                        pl.title(title)
                        divider = make_axes_locatable(pl.gca())
                        cax = divider.append_axes("right", "5%", pad="3%")
                        pl.colorbar(CS, cax=cax)
                    pl.tight_layout()
                    if save_pdf:
                        fname = "slice_equatorial_Ta{:g}".format(self.eq_params['taylor']) + ".pdf"
                        pl.savefig(fname, bbox_inches='tight', dpi=200)
                    if plot:
                        pl.show()
                    pl.clf()

    def apply_stencil(self, sol_vec):
        """Apply Galerkin stencil to recover physical fields"""
        
        for k,v in sol_vec.items():
            S = self.model.stencil(self.res, self.eq_params, self.eigs, self.bcs, k, False)
            sol_vec[k] = S*v

        return sol_vec

    def defaultWave(self, k):
        """Default conversion for wavenumber index"""

        return list(k)

def default_options():
    """Create dictionary of default options"""

    opts = dict()

    opts['fixed_shift'] = False     # Compute random shift only once
    opts['ellipse_radius'] = None   # Restrict eigenvalue search to be within radius
    opts['mode'] = 0                # Mode to track
    opts['root_tol'] = 1e-8         # Tolerance used in root finding algorighm
    opts['evp_tol'] = 1e-10         # Tolerance used in eigensolver

    opts['geometry'] = 'c1d'        # Geometry of the system

    opts['point'] = False           # Compute marginal curve point
    opts['point_k'] = 1           # Compute marginal curve point
    opts['curve'] = False           # Trace marginal curve
    opts['curve_points'] = [1]      # Wavenumbers for marginal curve
    opts['minimum'] = False         # Compute marginal curve minimum
    opts['minimum_int'] = False     # Only used integer wave numbers
    opts['solve'] = False           # Solve fixed point

    opts['plot_spy'] = False        # Plot matrix spy
    opts['plot_point'] = False      # Plot solution for marginal curve point
    opts['plot_curve'] = False      # Plot marginal curve 
    opts['viz_mode'] = 0        # Mode to visualize/save
    opts['show_spectra'] = False    # Plot solution spectra 
    opts['show_physical'] = False   # Plot solution in physical space 

    opts['write_mtx'] = False       # Write matrix to MatrixMarket files
    opts['save_curve'] = False      # Save marginal curve plot to pdf
    opts['save_spectra'] = False    # Save spectra plot to pdf
    opts['save_physical'] = False   # Save physical space plots to pdf
    opts['data_spectra'] = False    # Save spectral data to file(s)
    opts['data_physical'] = False   # Save physical data to file(s)

    return opts

def compute(gevp_opts, marginal_opts):
    """Compute eigen solutions requested by options"""

    # Make options consistent
    marginal_opts['minimum'] = marginal_opts['minimum'] and marginal_opts['curve']
    marginal_opts['plot_point'] = marginal_opts['plot_point'] and (marginal_opts['point'] or marginal_opts['minimum'])
    marginal_opts['plot_curve'] = marginal_opts['plot_curve'] and marginal_opts['curve']
    marginal_opts['show_spectra'] = marginal_opts['show_spectra'] and marginal_opts['solve']
    marginal_opts['show_physical'] = marginal_opts['show_physical'] and marginal_opts['solve']
    marginal_opts['save_spectra'] = marginal_opts['save_spectra'] and marginal_opts['solve']
    marginal_opts['save_physical'] = marginal_opts['save_physical'] and marginal_opts['solve']

    # Move options to GEVP
    gevp_opts['ellipse_radius'] = marginal_opts['ellipse_radius']
    gevp_opts['fixed_shift'] = marginal_opts['fixed_shift']

    if marginal_opts['point'] or marginal_opts['curve']:
        # Create marginal curve object
        curve = MarginalCurve(gevp_opts, mode = marginal_opts['mode'], rtol = marginal_opts['root_tol'], evp_tol = marginal_opts['evp_tol'])

    if marginal_opts['point'] and not marginal_opts['curve']:
        # Compute marginal curve at a single point
        kp = marginal_opts['point_k']
        Rac, evp_freq = curve.point(kp, guess = gevp_opts['eq_params']['rayleigh'])
        kc = kp

    if marginal_opts['curve']:
        # Trace marginal curve for a set of wave indexes
        ks = marginal_opts['curve_points']
        (data_k, data_Ra, data_freq) = curve.trace(ks, initial_guess = gevp_opts['eq_params']['rayleigh'])

        if marginal_opts['minimum']:
            # Compute minimum of marginal curve
            kc, Rac, fc = curve.minimum(data_k, data_Ra, only_int = marginal_opts['minimum_int'])

        if marginal_opts['plot_curve'] or marginal_opts['save_curve']:
            # Plot marginal curve and minimum
            if marginal_opts['minimum']:
                min_point = None
            else:
                min_point = (kc, Rac)
            curve.view(data_k, data_Ra, data_freq, minimum = min_point, plot = marginal_opts['plot_curve'], save_pdf = marginal_opts['save_curve'])

    if marginal_opts['plot_spy'] or marginal_opts['write_mtx'] or marginal_opts['solve']:
        if marginal_opts['plot_point']:
            Ra = Rac
            kp = kc
        else:
            kp = marginal_opts['point_k']
            Ra = gevp_opts['eq_params']['rayleigh']
        Print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(kp))
        gevp_opts['eigs'] = gevp_opts['wave'](kp)
        gevp_opts['tol'] = marginal_opts['evp_tol']
        gevp = GEVP(**gevp_opts)

    if marginal_opts['plot_spy'] or marginal_opts['write_mtx']:
        gevp.viewOperators(Ra, spy = marginal_opts['plot_spy'], write_mtx = marginal_opts['write_mtx'])

    if marginal_opts['solve']:
        gevp.solve(Ra, marginal_opts['viz_mode']+3, with_vectors = True)
        Print("Found eigenvalues:")
        Print(gevp.evp_lmb)

    if marginal_opts['show_spectra'] or marginal_opts['save_spectra'] or marginal_opts['data_spectra']:
        gevp.viewSpectra(marginal_opts['viz_mode'], plot = marginal_opts['show_spectra'], save_pdf = marginal_opts['save_spectra'])

    if marginal_opts['show_physical'] or marginal_opts['save_physical'] or marginal_opts['data_physical']:
        gevp.viewPhysical(marginal_opts['viz_mode'], marginal_opts['geometry'], plot = marginal_opts['show_physical'], save_pdf = marginal_opts['save_physical'])
