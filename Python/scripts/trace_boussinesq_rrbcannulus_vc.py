"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a cylindrical annulus model (velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rrbcannulus_vc as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRRBCAnnulusVC()
model.linearize = True
model.use_galerkin = False

## Set resolution, parameters, boundary conditions
#res = [14, 0, 14]
#eq_params = {'prandtl':7, 'rayleigh':3.440597e4, 'taylor':1e6, 'ro':2.0, 'rratio':0.75, 'scale3d':2.0}
#eigs = [0]
#eq_params = {'prandtl':7, 'rayleigh':3.56674e4, 'taylor':1e6, 'ro':2.0, 'rratio':0.75, 'scale3d':2.0}
#eq_params = {'prandtl':1, 'rayleigh':3.5e6, 'taylor':1e10, 'ro':1.0, 'rratio':0.35, 'scale3d':2.0}
#eq_params = {'prandtl':1, 'rayleigh':5e4, 'taylor':1e6, 'ro':1.0, 'rratio':0.35, 'scale3d':2.0}
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 2 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

# Create parameters
m = 7 
res = [64, 0, 64]
eq_params = {'taylor':1e6, 'prandtl':1, 'rayleigh':6e3, 'ro':1.0, 'rratio':0.35, 'scale3d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(0)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_point = True
marginal_curve = False
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (False and marginal_minimum)
marginal_show_point = (True and (marginal_point or marginal_minimum))
solve_gevp = (True or marginal_show_point)
show_spy = False
write_mtx = False
show_spectra = (True and solve_gevp) or marginal_show_point
show_physical = (True and solve_gevp) or marginal_show_point
viz_mode = 0

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts, rtol = 1e-8)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, evp_freq = curve.point(m, guess = eq_params['rayleigh'])

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ms = np.arange(max(0, m-5), m+6, 1)
    (data_m, data_Ra, data_freq) = curve.trace(ms, initial_guess = eq_params['rayleigh'])

    if marginal_minimum:
        # Compute minimum of marginal curve
        mc, Rac, fc = curve.minimum(data_m, data_Ra, only_int = True)

    if marginal_show_curve:
        # Plot marginal curve and minimum
        curve.view(data_m, data_Ra, data_freq, minimum = (mc, Rac), plot = True)

if show_spy or solve_gevp:
    if marginal_show_point:
        Ra = Rac
        m = mc
    else:
        Ra = eq_params['rayleigh']
    MarginalCurve.Print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(m))
    gevp_opts['eigs'] = wave(m)
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 1, with_vectors = True)
    MarginalCurve.Print("Found eigenvalues:")
    MarginalCurve.Print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode, naive = True)

if show_physical:
    gevp.viewPhysical(viz_mode, 'annulus', naive = True)
