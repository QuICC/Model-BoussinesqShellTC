"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcsphere as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCSphere()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
bc_vel = 1 # 0: NS, 1: SF
bc_temp = 0 # 0: FT 1: FF
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
#Ta = 1e8
#res = [32, 32, 0]
#Ta = 1e9
#res = [48, 48, 0]
#Ta = 1e10
#res = [64, 64, 0]
#Ta = 1e11
#res = [64, 64, 0]
Ta = 1e12
res = [96, 96, 0]
#Ta = 1e13
#res = [128, 128, 0]
#Ta = 1e14
#res = [256, 256, 0]
#Ta = 1e15
#res = [384, 384, 0]
#Ta = 1e16
#res = [512, 512, 0]
#Ta = 1e17
#res = [768, 768, 0]
#Ta = 1e18
#res = [1024, 1024, 0]

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
marginal_point = False
marginal_curve = True
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (False and marginal_minimum)
marginal_show_point = (True and (marginal_point or marginal_minimum))
solve_gevp = False or marginal_show_point
show_spy = False
write_mtx = False
show_spectra = (False and solve_gevp)
show_physical = (False and solve_gevp)
save_spectra = True
save_physical = True
viz_mode = 0

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts, rtol = 1e-8)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, evp_freq = curve.point(m, guess = eq_params['rayleigh'])
    mc = m

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ms = np.arange(max(1, m-5), m+6, 1)
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

if show_spectra or save_spectra:
    gevp.viewSpectra(viz_mode, plot = show_spectra, naive = True, save_pdf = save_spectra)

if show_physical or save_physical:
    gevp.viewPhysical(viz_mode, 'sphere', plot = show_physical, naive = True, save_pdf = save_physical)
