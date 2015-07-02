"""Script to run a marginal curve trace for the Boussinesq tilted F-plane 3DQG model"""

import numpy as np
import functools

import geomhdiscc.model.boussinesq_tilted_fplane3dqg_r as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTiltedFPlane3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [128, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'theta':0.0, 'scale1d':2.0}
kp = 1.3048
eq_params = {'prandtl':1, 'rayleigh':5.4780, 'theta':45.0, 'scale1d':2.0}
eq_params = {'prandtl':1, 'rayleigh':10.0, 'theta':45.0, 'scale1d':2.0}
eq_params = {'prandtl':1, 'rayleigh':1e-4, 'theta':45.0, 'scale1d':2.0}
kp = 1.1624
kp = 0.2199999999999
kp = 0.1

# Set wave number
phi = 90

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}

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
marginal_options['solve'] = True
marginal_options['point_k'] = kp
marginal_options['plot_point'] = True
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, kp-5), kp, kp+6)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
