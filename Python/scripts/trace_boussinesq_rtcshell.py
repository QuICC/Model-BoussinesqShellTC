"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcshell as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCShell()
model.linearize = True
model.use_galerkin = True

# Set resolution, parameters, boundary conditions
Rac = None
mc = None

# SF/SF, FT/FT, internal heating
#bc_vel = 1; bc_temp = 0; heating = 0; ro = 20./13.; rratio = 0.35
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
#Ta = 1e8; Rac = 31.957658460641; mc = 7
#res = [48, 48, 0]
#Ta = 1e9; Rac = 43.034758690274; mc = 10
#res = [48, 48, 0]
#Ta = 1e10; Rac = 59.975071459666; mc = 13
#res = [64, 64, 0]
#Ta = 1e11; Rac = 85.363356944817; mc = 20
#res = [64, 64, 0]
#Ta = 1e12; Rac = 122.69214718393; mc = 30
#res = [128, 128, 0]
#Ta = 1e13; Rac = 177.55422348123; mc = 44
#res = [192, 192, 0]
#Ta = 1e14; Rac = 258.13410447601; mc = 65
#res = [512, 256, 0]
#Ta = 1e15; Rac = 376.44742717745; mc = 95
#res = [512, 384, 0]
#Ta = 1e16
#res = [768, 512, 0]
#Ta = 1e17
#res = [768, 512, 0]
#Ta = 1e18
#res = [1024, 384, 0]

# NS/NS, FT/FT, internal heating
bc_vel = 0; bc_temp = 0; heating = 0; ro = 20./13.; rratio = 0.35
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
#Ta = 1e8; Rac = 31.534088376364; mc = 6 
#res = [48, 48, 0]
#Ta = 1e9; Rac = 42.219154540505; mc = 9
#res = [48, 48, 0]
#Ta = 1e10; Rac = 59.124359856967; mc = 13
#res = [96, 96, 0]
#Ta = 1e11; Rac = 84.487326687693; mc = 20
#res = [128, 128, 0]
#Ta = 1e12; Rac = 121.87739395205; mc = 30
#res = [192, 192, 0]
#Ta = 1e13; Rac = 176.79656879674; mc = 44
#res = [256, 256, 0]
#Ta = 1e14; Rac = 257.45628575047; mc = 65
res = [384, 256, 0]
Ta = 1e15; Rac = 375.86277729259; mc = 95
#res = [512, 384, 0]
#Ta = 1e16
#res = [768, 512, 0]
#Ta = 1e17
#res = [768, 768, 0]
#Ta = 1e18
#res = [1024, 1024, 0]

# Create parameters (rescaling to proper nondimensionalisation)
if mc is None:
    m = np.int(0.3029*Ta**(1./6.)) # Asymptotic prediction for minimum
else:
    m = mc
if Rac is None:
    Ra = (4.1173*Ta**(2./3.) + 17.7815*Ta**(0.5))*(1.0+rratio)/(2.0*ro**2*Ta**0.5) # Asymptotic prediction for critical Rayleigh number
else:
    Ra = Rac


res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
eq_params = {'taylor':Ta*(1.0-rratio)**4, 'prandtl':1, 'rayleigh':Ra, 'ro':ro, 'rratio':rratio, 'heating':heating}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(1)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-16
marginal_options['geometry'] = 'shell'
marginal_options['curve'] = False
marginal_options['minimum'] = False
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
