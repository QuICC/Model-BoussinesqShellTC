"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a plane layer (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np
import functools

import geomhdiscc.model.boussinesq_rrbcplane_vc as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqRRBCPlaneVC()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]

# FT
bc_temp = 0

# SF, FT,
bc_vel = 1
phi = 90
kp = 3.710
eq_params = {'prandtl':1, 'rayleigh':1676.12, 'taylor':0.0, 'heating':0, 'scale1d':2.0}
#kp = 129
#eq_params = {'prandtl':1, 'rayleigh':8.7050552e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 60
#eq_params = {'prandtl':1, 'rayleigh':4.04824521e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
## NS, FT,
#bc_vel = 0
#ky = 54
#eq_params = {'prandtl':1, 'rayleigh':3.455838e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55
#eq_params = {'prandtl':1, 'rayleigh':3.450289e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.4
#eq_params = {'prandtl':1, 'rayleigh':3.44979010e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.402
#eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 55.4023
#eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 56
#eq_params = {'prandtl':1, 'rayleigh':3.450893e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60
#eq_params = {'prandtl':1, 'rayleigh':3.515608e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 70
#eq_params = {'prandtl':1, 'rayleigh':4.134931e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 100
#eq_params = {'prandtl':1, 'rayleigh':1.094775e8, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#kp = 120
#eq_params = {'prandtl':1, 'rayleigh':7.7971364e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 121
#eq_params = {'prandtl':1, 'rayleigh':7.7892799e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 122
#eq_params = {'prandtl':1, 'rayleigh':7.7846430e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 123
#eq_params = {'prandtl':1, 'rayleigh':7.7832219e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 124
#eq_params = {'prandtl':1, 'rayleigh':7.7850143e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 128
#eq_params = {'prandtl':1, 'rayleigh':7.8420000e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 130
#eq_params = {'prandtl':1, 'rayleigh':7.8632272e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#kp = 131
#eq_params = {'prandtl':1, 'rayleigh':7.8875224e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generic Wave number function from single "index" (k perpendicular) and angle
def generic_wave(kp, phi):
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

# Wave number function from single "index" (k perpendicular)
wave = functools.partial(generic_wave, phi = phi)
eigs = wave(kp)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Setup computation, visualization and IO
marginal_point = False
marginal_curve = False
marginal_minimum = (True and marginal_curve)
marginal_show_curve = (True and marginal_minimum)
solve_gevp = True
show_spy = True
write_mtx = False
show_spectra = (True and solve_gevp)
show_physical = (True and solve_gevp)
viz_mode = 0

if marginal_point or marginal_curve:
    # Create marginal curve object
    curve = MarginalCurve.MarginalCurve(gevp_opts, rtol = 1e-8)

if marginal_point:
    # Compute marginal curve at a single point
    Rac, fc = curve.point(kp, initial_guess = 1e8)

if marginal_curve:
    # Trace marginal curve for a set of wave indexes
    ks = np.arange(121, 124, 1)
    (data_k, data_Ra, data_freq) = curve.trace(ks, guess = 7.8e8)

    if marginal_minimum:
        # Compute minimum of marginal curve
        kc, Rac, fc = curve.minimum(data_k, data_Ra)

    if marginal_show_curve:
        if marginal_minimum:
            minimum = (kc, Rac)
        else:
            minimum = None
        # Plot marginal curve and minimum (if available)
        curve.view(data_k, data_Ra, data_freq, minimum = minimum, plot = True)

if show_spy or solve_gevp:
    Ra = 0
    kp = 3.71
    print("Computing eigenvalues for Ra = " + str(Ra) + ", k = " + str(kp))
    gevp = MarginalCurve.GEVP(**gevp_opts)
    gevp.setEigs(kp)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 10, with_vectors = True)
    print("Found eigenvalues:")
    print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode)

if show_physical:
    gevp.viewPhysical(viz_mode)
