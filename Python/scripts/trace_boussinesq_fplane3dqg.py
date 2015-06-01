"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_fplane3dqg as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqFPlane3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'scale1d':2.0}
res = [64, 0, 0]
eq_params = {'prandtl':7, 'rayleigh':8.6957, 'scale1d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}
phi = 0

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(1.0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Create marginal curve object
curve = MarginalCurve.MarginalCurve(gevp_opts)

# Compute marginal curve at a single point
kp = 1.3048
Rac, evp_freq = curve.point(kp, guess = 5218.69901425)

# Trace marginal curve for a set of wave indexes
ks = np.arange(0.5, 10.5, 0.5)
(data_k, data_Ra, data_freq) = curve.trace(ks, initial_guess = 1)

## Compute minimum of marginal curve
kc, Rac, fc = curve.minimum(data_k, data_Ra)

# Plot marginal curve and minimum
curve.view(data_k, data_Ra, data_freq, minimum = (kc, Rac), plot = True)

# Setup computation, visualization and IO
solve_gevp = True
show_spy = False
show_spectra = (True and solve_gevp)
show_physical = (True and solve_gevp)
viz_mode = 0

if show_spy or solve_gevp:
    #kp = 1.5
    #Ra = 9.44899084506
    #kp = 3.0
    #Ra = 82.0966227126
    kp = 5.5
    Ra = 915.388767955
    kp = 8.5
    Ra = 5218.69901425
    print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(kp))
    gevp_opts['eigs'] = wave(kp)
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy:
    gevp.viewOperators(Ra, spy = True, write_mtx = True)

if solve_gevp:
    gevp.solve(Ra, 5, with_vectors = True)
    print("Found eigenvalues:")
    print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode)

if show_physical:
    gevp.viewPhysical(viz_mode)
