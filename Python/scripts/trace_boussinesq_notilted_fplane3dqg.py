"""Script to run a marginal curve trace for the Boussinesq tilted F-plane 3DQG model (nonorthogonal formulation)"""

import numpy as np
import functools

import geomhdiscc.model.boussinesq_notilted_fplane3dqg as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqNoTiltedFPlane3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]
kp = 1.3048529372573; Ra = 8.6956307662717
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'theta':0.0, 'scale1d':2.0}
#kp = 1.2898163457011; Ra = 8.3028296087578
#eq_params = {'prandtl':1, 'rayleigh':5.4780, 'theta':15.0, 'scale1d':2.0}
#kp = 1.2436965743243; Ra = 7.1780850234523
#eq_params = {'prandtl':1, 'rayleigh':4.0, 'theta':30.0, 'scale1d':2.0}
#kp = 1.1624855229293; Ra = 5.4779041217364
#eq_params = {'prandtl':1, 'rayleigh':4, 'theta':45.0, 'scale1d':2.0}

# Set wave number
phi = 0

bcs = {'bcType':model.SOLVER_HAS_BC, 'no_streamfunction':0, 'no_velocityz':0, 'temperature':0}

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

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['mode'] = 0
marginal_options['solve'] = False
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['point_k'] = kp
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = False
marginal_options['show_physical'] = False
marginal_options['curve_points'] = np.arange(max(1, kp-2), kp+3, 0.25)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
