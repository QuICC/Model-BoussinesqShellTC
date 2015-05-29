"""Module provides the functor required to trace a marginal curve"""

from __future__ import division
from __future__ import unicode_literals

import copy
import numpy as np
import geomhdiscc.linear_stability.solver_slepc as solver
import scipy.optimize as optimize

class MarginalCurve:
    """The marginal curve"""

    def __init__(self, gevp, mode = 0):
        """Initialize the marginal curve"""
        
        self.point = MarginalPoint(gevp, mode)

    def __call__(self, ks):
        """Trace the marginal curve"""
        
        nk = len(self.point.gevp.wave(ks[0]))
        data_k = np.zeros(len(ks))
        data_Ra = np.zeros(len(ks))
        data_freq = np.zeros(len(ks))
        Rac = 1
        for i, k in enumerate(ks):
            self.point.gevp.setEigs(k)
            Rac, freq = self.point(Rac)
            data_k[i] = k
            data_Ra[i] = Rac
            data_freq[i] = freq

        return (data_k, data_Ra, data_freq)

    def minimum(self, ks, Ras):
        """Compute the critical value (minimum of curve)"""
        
        # Get bounds for minimum
        imin = np.argmin(Ras)
        a = ks[imin-1]
        b = ks[imin]
        c = ks[imin+1]

        # Find minimum
        opt = optimize.minimize_scalar(self.critical, [a, b, c])

        return (opt.x, opt.fun, self.point.frequency())

    def critical(self, k):

        self.point.gevp.setEigs(k)
        return self.point(1)[0]

class MarginalPoint:
    """A point on the marginal curve"""

    def __init__(self, gevp, mode = 0):
        """Initialize the marginal point"""

        self.gevp = gevp
        self.mode = mode
        self.Rac = None
        self.fc = None

    def __call__(self, guess, rtol = 1e-7):
        """Compute marginal point in given interval [a, b]"""

        a = 1
        b = 100
        Rac = optimize.brentq(self.tracker, a, b, rtol = rtol)
        self.Rac = Rac
        self.fc = self.frequency()

        return (self.Rac, self.fc)

    def tracker(self, Ra):
        """Track the growth rate to find critical value"""

        self.gevp.solve(Ra, self.mode+1)

        return self.gevp.evp_lmb[self.mode].real

    def eigenvector(self):
        """Compute the eigenvectors"""

        self.gevp.solve(self.Rac, self.mode+1, with_vectors = True)

        return self.gevp.evp_vec[:,self.mode]

    def eigenvalue(self):
        """Compute the eigenvalues"""

        self.gevp.solve(self.Rac, self.mode+1)

        return self.growthRate() + 1j*self.frequency()

    def growthRate(self):
        """Get the growth rate"""

        self.gevp.solve(self.Rac, self.mode+1)

        return self.gevp.evp_lmb[self.mode].real

    def frequency(self):
        """Compute the growth rate"""

        self.gevp.solve(self.Rac, self.mode+1)

        return self.gevp.evp_lmb[self.mode].imag

class GEVP:
    """Class to represent a marginal point on a marginal curve"""
    
    def __init__(self, model, res, eq_params, eigs, bcs, fields, wave = None):
        """Initialize the marginal point variables"""

        self.model = model
        self.res = res
        self.eq_params = eq_params.copy()
        self.eigs = eigs
        self.bcs = bcs.copy()
        self.no_bcs = bcs.copy()
        self.no_bcs['bcType'] = model.SOLVER_NO_TAU
        self.fields = fields
        self.evp_lmb = None
        self.evp_vec = None
        if wave is None:
            self.wave = self.defaultWave
        else:
            self.wave = wave

    def buildMatrices(self, Ra):
        """Build A and B matrices for the given rayleigh number"""

        self.eq_params['rayleigh'] = Ra
        A = self.model.implicit_linear(self.res, self.eq_params, self.eigs, self.bcs, self.fields)
        B = self.model.time(self.res, self.eq_params, self.eigs, self.no_bcs, self.fields)

        return (A,B)
   
    def setEigs(self, k):
       """Set the wave number"""

       self.eigs = self.wave(k)

    def solve(self, Ra, nev, with_vectors = False):
        """Solve the GEVP and store eigenvalues"""

        if Ra != self.eq_params['rayleigh'] or (with_vectors and self.evp_vec is None):
            A, B = self.buildMatrices(Ra)
            if with_vectors:
                self.evp_lmb, self.evp_vec = solver.eigenpairs(A, B, nev)
            else:
                self.evp_lmb = solver.eigenvalues(A, B, nev)
                self.evp_vec = None

    def defaultWave(self, k):
        return list(k)
