"""Script to run a marginal curve trace for the Boussinesq thermal convection in a sphere (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_tcsphere_std as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqTCSphereStd()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]
l = 7; eq_params = {'prandtl':1, 'rayleigh':1.645e4}
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
marginal_point = True
marginal_curve = False
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (False and marginal_minimum)
marginal_show_point = (True and (marginal_point or marginal_minimum))
solve_gevp = (True or marginal_show_point)
show_spy = True
write_mtx = False
show_spectra = (True and solve_gevp) or marginal_show_point
show_physical = (False and solve_gevp) or marginal_show_point
viz_mode = 0

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts, rtol = 1e-8)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, evp_freq = curve.point(l, guess = eq_params['rayleigh'])
    lc = l

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ms = np.arange(max(0, l-5), l+6, 1)
    (data_l, data_Ra, data_freq) = curve.trace(ms, initial_guess = eq_params['rayleigh'])

    if marginal_minimum:
        # Compute minimum of marginal curve
        lc, Rac, fc = curve.minimum(data_l, data_Ra, only_int = True)

    if marginal_show_curve:
        # Plot marginal curve and minimum
        curve.view(data_l, data_Ra, data_freq, minimum = (lc, Rac), plot = True)

if show_spy or solve_gevp:
    if marginal_show_point:
        Ra = Rac
        l = lc
    else:
        Ra = eq_params['rayleigh']
    MarginalCurve.Print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(l))
    gevp_opts['eigs'] = wave(l)
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 1, with_vectors = True)
    MarginalCurve.Print("Found eigenvalues:")
    MarginalCurve.Print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode)

if show_physical:
    gevp.viewPhysical(viz_mode, 'b1d')
