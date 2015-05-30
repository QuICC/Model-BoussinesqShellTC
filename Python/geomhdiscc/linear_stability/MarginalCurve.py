"""Module provides the functor required to trace a marginal curve"""

from __future__ import division
from __future__ import unicode_literals

import copy
import numpy as np
import geomhdiscc.linear_stability.solver_slepc as solver
import scipy.optimize as optimize

class NewtonNoneError(Exception):
    pass

class NewtonDoneError(Exception):
    pass

class MarginalCurve:
    """The marginal curve"""

    def __init__(self, gevp_opts, mode = 0):
        """Initialize the marginal curve"""
        
        self.point = MarginalPoint(gevp_opts, mode)

    def trace(self, ks):
        """Trace the marginal curve"""
        
        data_k = np.zeros(len(ks))
        data_Ra = np.zeros(len(ks))
        data_freq = np.zeros(len(ks))
        Rac = None
        for i, k in enumerate(ks):
            if i > 3:
                Rac = p(k)
            Rac, freq = self.point(k, guess = Rac)
            data_k[i] = k
            data_Ra[i] = Rac
            data_freq[i] = freq
            if i > 2:
                lb = max(0, i-3)
                z = np.polyfit(data_k[lb:i+1], data_Ra[lb:i+1], min(i,3))
                p = np.poly1d(z)

        return (data_k, data_Ra, data_freq)

    def minimum(self, ks, Ras):
        """Compute the critical value (minimum of curve)"""
        
        # Get bounds for minimum
        imin = np.argmin(Ras)
        a = ks[imin-1]
        b = ks[imin]
        c = ks[imin+1]

        # Find minimum
        opt = optimize.minimize_scalar(self._tracker, [a, b, c])

        return (opt.x, opt.fun, self.point.frequency())

    def _tracker(self, k):

        return self.point(k)[0]

class MarginalPoint:
    """A point on the marginal curve"""

    def __init__(self, gevp_opts, mode = 0):
        """Initialize the marginal point"""

        self._gevp = GEVP(**gevp_opts)
        self.mode = mode
        self.Rac = None
        self.fc = None
        self.a = None
        self.b = None

    def __call__(self, k, guess = None, rtol = 1e-7):
        """Compute marginal point in given interval [a, b]"""

        if guess is None:
            Ra = 1
            factor = 2
        else:
            Ra = guess
            factor = 1.1

        print("Computing marginal curve for k = " + str(k))
        print("\t- initial guess Ra = " + str(Ra))

        # Set wavenumter
        self._gevp.setEigs(k)

        # Search interval (newton)
        print("\t- finding interval")
        for i in range(0,16):
            self.findInterval(Ra)
            if self.a != None:
                break
            else:
                Ra = factor*Ra

        if self.a is None:
            self.Rac = float('NaN')
            self.rc = float('NaN')
            self.a = None
            self.b = None
        else:
            # Refine zero with Brent's method
            print("\t- Brent's method: [" + str(self.a) + ', ' + str(self.b) + "]")
            Rac = optimize.brentq(self._tracker, self.a, self.b, rtol = rtol)
            self.a = None
            self.b = None
            self.Rac = Rac
            self.fc = self.frequency()
        print("\t- Ra = " + str(self.Rac) + ", freq = " + str(self.fc))

        return (self.Rac, self.fc)

    def findInterval(self, guess):
        """Find the sign changing interval"""

        try:
            Rac = optimize.newton(self._signTracker, guess)
        except (NewtonDoneError, NewtonNoneError):
            pass
        except:
            raise

    def eigenvector(self):
        """Compute the eigenvectors"""

        self._gevp.solve(self.Rac, self.mode+1, with_vectors = True)

        return self._gevp.evp_vec[:,self.mode]

    def eigenvalue(self):
        """Compute the eigenvalues"""

        self._gevp.solve(self.Rac, self.mode+1)

        return self._gevp.evp_lmb[self.mode]

    def growthRate(self):
        """Get the growth rate"""

        self._gevp.solve(self.Rac, self.mode+1)

        return self._gevp.evp_lmb[self.mode].real

    def frequency(self):
        """Compute the growth rate"""

        self._gevp.solve(self.Rac, self.mode+1)

        return self._gevp.evp_lmb[self.mode].imag

    def _tracker(self, Ra):
        """Track the growth rate to find critical value"""

        self._gevp.solve(Ra, self.mode+1)

        return self._gevp.evp_lmb[self.mode].real

    def _signTracker(self, Ra):
        """Track the growth rate to find critical value"""

        print(Ra)
        # Get change in rayleigh number
        dRa = Ra - self._gevp.eq_params['rayleigh']
        # Solve
        self._gevp.solve(Ra, self.mode+1)
        print(self._gevp.evp_lmb)

        if self._gevp.evp_lmb is None:
            raise NewtonNoneError

        growth = self._gevp.evp_lmb[self.mode].real
        # Check for change of sign (require small change in rayleigh and small growth rate)
        if abs(dRa/Ra) < 0.1 and abs(self._gevp.evp_lmb[self.mode].real) < 1:
            self._gevp.solve(Ra+dRa, self.mode+1)
            if growth*self._gevp.evp_lmb[self.mode].real < 0:
                self.a = min(Ra, Ra+dRa)
                self.b = max(Ra, Ra+dRa)
                raise NewtonDoneError

        return growth

class GEVP:
    """Class to represent a marginal point on a marginal curve"""
    
    def __init__(self, model, res, eq_params, eigs, bcs, wave = None):
        """Initialize the marginal point variables"""

        self.model = copy.copy(model)
        self.res = res
        self.eq_params = eq_params.copy()
        self.eigs = eigs
        self.bcs = bcs.copy()
        self.no_bcs = bcs.copy()
        self.no_bcs['bcType'] = model.SOLVER_NO_TAU
        self.fields = self.model.stability_fields()
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
        """Default conversion for wavenumber index"""

        return list(k)
