"""Script to run a marginal curve trace for the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_tcshell_std as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTCShellStd()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]
#l = 1; eq_params = {'prandtl':1, 'rayleigh':1.645e4, 'ro':1, 'rratio':0.2}
l = 1; eq_params = {'prandtl':1, 'rayleigh':3.0468e4, 'ro':1, 'rratio':0.3, 'heating':0} #NS/NS
#l = 1; eq_params = {'prandtl':1, 'rayleigh':8.4791e3, 'ro':1, 'rratio':0.3} #SF/SF
#l = 1; eq_params = {'prandtl':1, 'rayleigh':4.183202e4, 'ro':1, 'rratio':0.5} #SF/SF
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

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
marginal_options['geometry'] = 's1d'
marginal_options['curve'] = True
marginal_options['minimum'] = True
marginal_options['plot_curve'] = True
marginal_options['solve'] = True
marginal_options['minimum_int'] = True
marginal_options['point_k'] = l
marginal_options['plot_point'] = True
marginal_options['show_spectra'] = True
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, l-2), l+3, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)
