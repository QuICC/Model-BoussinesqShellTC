"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcshell as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCShell()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
#eq_params = {'taylor':1e10, 'prandtl':1, 'rayleigh':2.073175e7, 'ro':1, 'rratio':0.35} # m = 13, NS/NS
#eq_params = {'taylor':1e10, 'prandtl':1, 'rayleigh':2.103005e7, 'ro':1, 'rratio':0.35} # m = 13, SF/SF
#eq_params = {'taylor':1e12, 'prandtl':1, 'rayleigh':4.2726e8, 'ro':1, 'rratio':0.35} # m = 30, NS/NS
#eq_params = {'taylor':1e12, 'prandtl':1, 'rayleigh':4.302e8, 'ro':1, 'rratio':0.35} # m = 30, SF/SF
#eq_params = {'taylor':1e10, 'prandtl':1, 'rayleigh':3.2e6, 'ro':20.0/13.0, 'rratio':0.35} # m = 4, NS/NS
#eq_params = {'taylor':1e10, 'prandtl':1, 'rayleigh':2.07310821e7, 'ro':1.0, 'rratio':0.35, 'heating':0} # m = 13, NS/NS, internal heating
#eq_params = {'taylor':1e10, 'prandtl':1, 'rayleigh':3.24308e6, 'ro':1.0, 'rratio':0.35, 'heating':1} # m = 9, NS/NS, differential heating
#eq_params = {'taylor':4e6, 'prandtl':1, 'rayleigh':4.75992e4, 'ro':1.0, 'rratio':0.35, 'heating':1} # m = 9, NS/NS, differential heating
#eq_params = {'taylor':4e6*0.65**4, 'prandtl':1, 'rayleigh':28.7296, 'ro':20./13., 'rratio':0.35, 'heating':1} # m = 9, NS/NS, differential heating
#eq_params = {'taylor':1e10*(1.0-0.35)**4, 'prandtl':1, 'rayleigh':3.24308e6/((20./13.)**2*0.35*1e10**0.5), 'ro':20./13., 'rratio':0.35, 'heating':1} # m = 9, NS/NS, differential heating
# Internal heating, velocity bc: NS/NS, temperature bc: FT/FT
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
ro = 20./13.
rratio = 0.35
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
#res = [192, 192, 0]
#Ta = 1e15
#res = [256, 256, 0]
#Ta = 1e16
#res = [384, 384, 0]
#Ta = 1e17
#res = [384, 384, 0]
#Ta = 1e18
#res = [512, 512, 0]

# Create parameters (rescaling to proper nondimensionalisation)
m = np.int(0.3029*Ta**(1./6.)) # Asymptotic prediction for minimum
res = [res[0], res[1]+m, 0] # Extend harmonic degree by harmonic order (fixed number of modes)
Ra_th = 4.1173*Ta**(2./3.) + 17.7815*Ta**(0.5) # Asymptotic prediction for critical Rayleigh number
eq_params = {'taylor':Ta*(1.0-rratio)**4, 'prandtl':1, 'rayleigh':Ra_th*(1.0+rratio)/(2.0*ro**2*Ta**0.5), 'ro':ro, 'rratio':rratio, 'heating':0}
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
solve_gevp = True
show_spy = False
write_mtx = False
show_spectra = (True and solve_gevp)
show_physical = (True and solve_gevp)
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
    gevp.viewPhysical(viz_mode, 'shell', naive = True)
