"""Script to run a marginal curve trace for the Boussinesq inertial wave in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_iwshell as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqIWShell()
model.linearize = True
model.use_galerkin = True

# Set resolution, parameters, boundary conditions
mc = None

# SF/SF
#bc_vel = 1; ro = 20./13.; rratio = 0.35
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]

# NS/NS
bc_vel = 0; ro = 20./13.; rratio = 0.35
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
Ta = 1e12; mc = 0
res = [128, 128, 0]

# Create parameters (rescaling to proper nondimensionalisation)
if mc is None:
    m = np.int(0.3029*Ta**(1./6.)) # Asymptotic prediction for minimum
else:
    m = mc

res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
eq_params = {'taylor':Ta*(1.0-rratio)**4, 'prandtl':1, 'ro':ro, 'rratio':rratio, 'rayleigh':0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel}

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
marginal_options['plot_spy'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['viz_mode'] = 3
marginal_options['curve_points'] = np.arange(max(0, m-2), m+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
