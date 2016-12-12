"""Script to run a marginal curve trace for the Boussinesq inertial waves in a cylinder model with Worland expansion (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_iwcylinder as mod
import geomhdiscc.linear_stability.marginal_curve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqIWCylinder()
model.linearize = True
model.use_galerkin = False

# Set boundary conditions
bc_vel = 0 # 0: NS/NS

# Create parameters
m = 3
res = [64, 0, 64]
Ta = (1e-3)**(-2)
Gamma = 1.0
eq_params = {'taylor':Ta, 'gamma':Gamma, 'rayleigh':0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_options = MarginalCurve.default_options()
marginal_options['evp_tol'] = 1e-12
marginal_options['geometry'] = 'cylinder_worland'
marginal_options['ellipse_radius'] = 1e3
marginal_options['target'] = -104
marginal_options['curve'] = False
marginal_options['minimum'] = False
marginal_options['minimum_int'] = True
marginal_options['plot_curve'] = False
marginal_options['solve'] = True
marginal_options['solve_nev'] = 3
marginal_options['point_k'] = m
marginal_options['plot_point'] = False
marginal_options['plot_spy'] = True
marginal_options['write_mtx'] = True
marginal_options['show_spectra'] = True
marginal_options['viz_mode'] = -1
marginal_options['show_physical'] = True
marginal_options['curve_points'] = np.arange(max(0, m-0), m+1, 1)

# Compute 
MarginalCurve.compute(gevp_opts, marginal_options)