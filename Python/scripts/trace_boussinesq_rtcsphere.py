"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcsphere as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCSphere()
model.linearize = True
model.use_galerkin = True

# Set resolution, parameters, boundary conditions
bc_vel = 1 # 0: NS, 1: SF
bc_temp = 0 # 0: FT 1: FF
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
Ta = 1e8
res = [32, 32, 0]
#Ta = 1e9
#res = [48, 48, 0]
#Ta = 1e10
#res = [64, 64, 0]
#Ta = 1e11
#res = [96, 96, 0]
#Ta = 1e12
#res = [128, 128, 0]
#Ta = 1e13
#res = [192, 192, 0]
#Ta = 1e14
#res = [256, 256, 0]
#Ta = 1e15
#res = [512, 512, 0]
#Ta = 1e16
#res = [512, 768, 0]
#Ta = 1e17
#res = [512, 1024, 0]
#Ta = 1e18
#res = [784, 1536, 0]
#Ta = 1e19
#res = [512, 512, 0]

# Create parameters (rescaling to proper nondimensionalisation)
m = np.int(0.3029*Ta**(1./6.)) # Asymptotic prediction for minimum
res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
Ra_th = 4.1173*Ta**(2./3.) + 17.7815*Ta**(0.5) # Asymptotic prediction for critical Rayleigh number
eq_params = {'taylor':Ta, 'prandtl':1, 'rayleigh':Ra_th}
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
marginal_options['geometry'] = 'sphere'
marginal_options['curve'] = False
marginal_options['minimum'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['minimum_int'] = True
marginal_options['point_k'] = m
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
