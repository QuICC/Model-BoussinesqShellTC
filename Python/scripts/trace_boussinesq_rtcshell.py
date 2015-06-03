"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcshell as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRTCShell()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
l = 0
m = 300
res = [512, m+512, 0]
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
eq_params = {'taylor':1e10*(1.0-0.35)**4, 'prandtl':1, 'rayleigh':2.07310821e7*(1.0+0.35)/(2.0*(20./13.)**2*1e10**0.5), 'ro':20./13., 'rratio':0.35, 'heating':0} # m = 13, NS/NS, internal heating
eq_params = {'taylor':1e12*(1.0-0.35)**4, 'prandtl':1, 'rayleigh':2.07310821e7*(1.0+0.35)/(2.0*(20./13.)**2*1e12**0.5), 'ro':20./13., 'rratio':0.35, 'heating':0} # m = 13, NS/NS, internal heating
eq_params = {'taylor':1e14*(1.0-0.35)**4, 'prandtl':1, 'rayleigh':2.07310821e7*(1.0+0.35)/(2.0*(20./13.)**2*1e14**0.5), 'ro':20./13., 'rratio':0.35, 'heating':0} # m = 13, NS/NS, internal heating
eq_params = {'taylor':1e18*(1.0-0.35)**4, 'prandtl':1, 'rayleigh':2.07310821e7*(1.0+0.35)/(2.0*(20./13.)**2*1e18**0.5), 'ro':20./13., 'rratio':0.35, 'heating':0} # m = 13, NS/NS, internal heating
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(13)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_point = False
marginal_curve = True
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (False and marginal_minimum)
solve_gevp = False
show_spy = False
write_mtx = False
show_spectra = (True and solve_gevp)
show_physical = (False and solve_gevp)
viz_mode = 0

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, evp_freq = curve.point(m, guess = 600)

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ms = np.arange(300, 310, 1)
    (data_m, data_Ra, data_freq) = curve.trace(ms, initial_guess = 650)

    if marginal_minimum:
        # Compute minimum of marginal curve
        mc, Rac, fc = curve.minimum(data_m, data_Ra, only_int = True)

    if marginal_show_curve:
        # Plot marginal curve and minimum
        curve.view(data_m, data_Ra, data_freq, minimum = (mc, Rac), plot = True)

if show_spy or solve_gevp:
    Ra = 1000
    print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(m))
    gevp_opts['eigs'] = wave(m)
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 5, with_vectors = True)
    print("Found eigenvalues:")
    print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode, naive = True)

if show_physical:
    gevp.viewPhysical(viz_mode)
