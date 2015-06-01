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
m = 13
res = [128, 128, 0]
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
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Wave number function from single "index" (k perpendicular)
def wave(m):
    return [float(m)]

eigs = wave(13)

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Create marginal curve object
curve = MarginalCurve.MarginalCurve(gevp_opts)

# Compute marginal curve at a single point
#m = 13
#Rac, evp_freq = curve.point(m, guess = 2e7)
#print((kp, Rac, evp_freq))

## Trace marginal curve for a set of wave indexes
#ms = np.arange(1, 11, 1)
#(data_m, data_Ra, data_freq) = curve.trace(ms)
#
## Compute minimum of marginal curve
#mc, Rac, fc = curve.minimum(data_m, data_Ra)
#print((mc, Rac, fc))
#
## Plot marginal curve and minimum
#import matplotlib.pylab as pl
#pl.plot(data_m, data_Ra, 'b', data_m, data_freq, 'g', mc, Rac, 'r+', markersize=14)
#pl.show()

# Setup comutation, visualization and IO
solve_gevp = True
show_spy = False
write_mtx = False
show_spectra = (True and solve_gevp)
show_physical = (True and solve_gevp)
viz_mode = 0

if show_spy or solve_gevp:
    m = 13
    Ra = 2e7
    gevp_opts['eigs'] = wave(m)
    gevp = MarginalCurve.GEVP(**gevp_opts)

if show_spy or write_mtx:
    gevp.viewOperators(Ra, spy = show_spy, write_mtx = write_mtx)

if solve_gevp:
    gevp.solve(Ra, 5, with_vectors = True)
    print(gevp.evp_lmb)

if show_spectra:
    gevp.viewSpectra(viz_mode)

if show_physical:
    gevp.viewPhysical(viz_mode)
#
#
#
#
## Setup visualization and IO
#show_spy = False
#write_mtx = True
#solve_evp = True
#show_solution = (True and solve_evp)
#
#if show_spy or show_solution:
#    import matplotlib.pylab as pl
#
#if show_solution:
#    import geomhdiscc.transform.shell as transf
#
## Show the "spy" of the two matrices
#if show_spy:
#    pl.spy(A, markersize=3, marker = '.', markeredgecolor = 'b')
#    pl.tick_params(axis='x', labelsize=30)
#    pl.tick_params(axis='y', labelsize=30)
#    pl.show()
#    pl.spy(B, markersize=3, marker = '.', markeredgecolor = 'b')
#    pl.tick_params(axis='x', labelsize=30)
#    pl.tick_params(axis='y', labelsize=30)
#    pl.show()
#
## Export the two matrices to matrix market format
#if write_mtx:
#    import scipy.io as io
#    io.mmwrite("matrix_A.mtx", A)
#    io.mmwrite("matrix_B.mtx", B)
#    print("Wrote matrices to MTX files")
#
## Solve EVP with sptarn
#if solve_evp:
#    import geomhdiscc.linear_stability.solver as solver
#    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, 1e2)
#    print(evp_lmb)
#
#if show_solution:
#    viz_mode = 0
#    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
#    # Get solution vectors
#    nL = res[1]-m
#    if model.use_galerkin == True:
#        nTor = 2
#        nPol = 4
#        nT = 2
#    else:
#        nTor = 0
#        nPol = 0
#        nT = 0
#    sol_tor = evp_vec[0:(res[0]-nTor)*nL,viz_mode]
#    sol_pol = evp_vec[(res[0]-nTor)*nL:(2*res[0]-nTor-nPol)*nL,viz_mode]
#    sol_t = evp_vec[(2*res[0]-nTor-nTor)*nL:(3*res[0]-nTor-nPol-nT)*nL,viz_mode]
#    
#    # Create spectrum plots
#    pl.subplot(1,3,1)
#    pl.semilogy(np.abs(sol_tor))
#    pl.title("Toroidal")
#    pl.subplot(1,3,2)
#    pl.semilogy(np.abs(sol_pol))
#    pl.title("Poloidal")
#    pl.subplot(1,3,3)
#    pl.semilogy(np.abs(sol_t))
#    pl.title("T")
#    pl.show()
#    pl.close("all")
