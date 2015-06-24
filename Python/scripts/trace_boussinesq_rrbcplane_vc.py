"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np
import functools

import geomhdiscc.model.boussinesq_rrbcplane_vc as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRRBCPlaneVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
Rac = None
kc = None

# SF, FT,
#bc_vel = 1
#heating = 0
#bc_temp = 0
#phi = 0

#phi = 45
#kp = 3.710
#kp = 129
#kp = 60

# NS, FT,
bc_vel = 0
heating = 0
bc_temp = 0
phi = 0

res = [64, 0, 0]
Ta = 1e8;
kp = 20

#eq_params = {'prandtl':1, 'rayleigh':2.1544e6, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 54
#eq_params = {'prandtl':1, 'rayleigh':3.455838e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55
#eq_params = {'prandtl':1, 'rayleigh':3.450289e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.4
#eq_params = {'prandtl':1, 'rayleigh':3.44979010e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.402
#eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.4023
#eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 56
#eq_params = {'prandtl':1, 'rayleigh':3.450893e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60
#eq_params = {'prandtl':1, 'rayleigh':3.515608e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 70
#eq_params = {'prandtl':1, 'rayleigh':4.134931e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 100
#eq_params = {'prandtl':1, 'rayleigh':1.094775e8, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 120
#eq_params = {'prandtl':1, 'rayleigh':7.7971364e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 121
#eq_params = {'prandtl':1, 'rayleigh':7.7892799e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 122
#eq_params = {'prandtl':1, 'rayleigh':7.7846430e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 123
#eq_params = {'prandtl':1, 'rayleigh':7.7832219e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 124
#eq_params = {'prandtl':1, 'rayleigh':7.7850143e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 128
#eq_params = {'prandtl':1, 'rayleigh':7.8420000e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 130
#eq_params = {'prandtl':1, 'rayleigh':7.8632272e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 131
#eq_params = {'prandtl':1, 'rayleigh':7.8875224e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}

# Create parameters (rescaling to proper nondimensionalisation)
if kc is None:
    kp = Ta**(1./6.) # Asymptotic prediction for minimum
else:
    kp = kc
if Rac is None:
    Ra = Ta**(2./3.) # Asymptotic prediction for critical Rayleigh number
else:
    Ra = Rac

eq_params = {'taylor':Ta, 'prandtl':1, 'rayleigh':Ra, 'heating':heating, 'scale1d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(kp)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['ellipse_radius'] = 1e5
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = kp
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, kp-5), kp+6, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
