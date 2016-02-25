"""Script to run a marginal curve trace for the Boussinesq thermal convection in a sphere (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_tcsphere_std as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTCSphereStd()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#res = [64, 0, 0]
#bc_vel = 0 # NS/NS
#bc_temp = 0 # FT/FT
#l = 1; eq_params = {'prandtl':1, 'rayleigh':1.645e4}

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]
bc_vel = 1 # SF/SF
bc_temp = 0 # FT/FT
l = 1; eq_params = {'prandtl':1, 'rayleigh':1.645e4}

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(l):
    return [float(l)]

eigs = wave(l)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-16
marginal_options['geometry'] = 'b1d'
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = False
marginal_options['point_k'] = l
marginal_options['plot_point'] = True
marginal_options['plot_spy'] = False
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(1, l-10), l+11, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
