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
bc_vel = 1; bc_temp = 0; heating = 0; ro = 20./13.; rratio = 0.35
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
res = [512, 256, 0]
Ta = 1e15; Rac = 376.44742717745; mc = 95
#res = [512, 384, 0]
#Ta = 1e16
#res = [768, 512, 0]
#Ta = 1e17
#res = [768, 512, 0]
#Ta = 1e18
#res = [1024, 384, 0]

# NS/NS, FT/FT, internal heating
#bc_vel = 0; bc_temp = 0; heating = 0; ro = 20./13.; rratio = 0.35
#Ta = 1e6
#res = [32, 32, 0]
#Ta = 1e7
#res = [32, 32, 0]
#Ta = 1e8
#res = [48, 48, 0]
#Ta = 1e9
#res = [48, 48, 0]
#Ta = 1e10
#res = [64, 64, 0]
#Ta = 1e11
#res = [64, 64, 0]
#Ta = 1e12
#res = [128, 128, 0]
#Ta = 1e13
#res = [192, 192, 0]
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
marginal_point = False
marginal_curve = False
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (False and marginal_minimum)
marginal_show_point = (False and (marginal_point or marginal_minimum))
solve_gevp = (True or marginal_show_point)
show_spy = False
write_mtx = False
show_spectra = (True and solve_gevp) or marginal_show_point
show_physical = (False and solve_gevp) or marginal_show_point
save_spectra = False
save_physical = False
viz_mode = 0

evp_tol = 1e-30

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts, rtol = 1e-8, evp_tol = evp_tol)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, evp_freq = curve.point(m, guess = eq_params['rayleigh'])
    mc = m

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ms = np.arange(max(0, m-3), m+4, 1)
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
    gevp_opts['tol'] = evp_tol
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 1, with_vectors = True)
    MarginalCurve.Print("Found eigenvalues:")
    MarginalCurve.Print(gevp.evp_lmb)

if solve_gevp:
    gevp.viewSpectra(viz_mode, plot = show_spectra, save_pdf = save_spectra)
    gevp.viewPhysical(viz_mode, 'shell', plot = show_physical, save_pdf = save_physical)
