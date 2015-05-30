"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_fplane3dqg as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqFPlane3DQG()
model.linearize = True
model.use_galerkin = False

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'scale1d':2.0}
res = [64, 0, 0]
eq_params = {'prandtl':7, 'rayleigh':8.6957, 'scale1d':2.0}
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}

# Wave number function from single "index" (k perpendicular)
def wave(kp):
    phi = 0
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

eigs = [1, 0]

# Collect GEVP setup parameters into single dictionary
gevp_opts = {'model':model, 'res':res, 'eq_params':eq_params, 'eigs':eigs, 'bcs':bcs, 'wave':wave}

# Create marginal curve object
curve = MarginalCurve.MarginalCurve(gevp_opts)

# Compute marginal curve at a single point
#kp = 1.3048
#kp = 25.0
#Rac, evp_freq = curve.point(kp)
#print((kp, Rac, evp_freq))
#
## Trace marginal curve for a set of wave indexes
#ks = np.arange(0.5, 100.5, 0.5)
#(data_k, data_Ra, data_freq) = curve.trace(ks)
#
## Compute minimum of marginal curve
#kc, Rac, fc = curve.minimum(data_k, data_Ra)
#print((kc, Rac, fc))
#
# Plot marginal curve and minimum
#import matplotlib.pylab as pl
#pl.plot(data_k, data_Ra, 'b', data_k, data_freq, 'g', kc, Rac, 'r+', markersize=14)
#pl.show()

# Compute and visualize some mode
kp = 7
Ra = 2402
gevp_opts['eigs'] = wave(kp)
gevp = MarginalCurve.GEVP(**gevp_opts)

# Setup visualization and IO
show_spy = False
write_mtx = False
solve_gevp = True
show_solution = (True and solve_gevp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

if show_spy or write_mtx:
    A, B = gevp.buildMatrices(Ra)

# Show the "spy" of the two matrices
if show_spy:
    pl.spy(A, markersize=5, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()
    pl.spy(B, markersize=5, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve GEVP
if solve_gevp:
    gevp.solve(Ra, 5, with_vectors = True)
    evp_lmb = gevp.evp_lmb
    evp_vec = gevp.evp_vec
    print(evp_lmb)

if show_solution:
    viz_mode = 0
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    sol_s = evp_vec[0:res[0],viz_mode]
    sol_w = evp_vec[res[0]:2*res[0],viz_mode]
    sol_t = evp_vec[2*res[0]:3*res[0],viz_mode]
    # Create spectrum plots
    pl.subplot(1,3,1)
    pl.semilogy(abs(sol_s))
    pl.title('s')
    pl.subplot(1,3,2)
    pl.semilogy(abs(sol_w))
    pl.title('w')
    pl.subplot(1,3,3)
    pl.semilogy(abs(sol_t))
    pl.title('T')
    pl.show()

    # Compute physical space values
    grid_x = transf.grid(res[0])
    phys_s = transf.tophys(sol_s.real)
    phys_w = transf.tophys(sol_w.real)
    phys_t = transf.tophys(sol_t.real)

    # Show physical plot
    pl.subplot(1,3,1)
    pl.plot(grid_x, phys_s)
    pl.title('s')
    pl.subplot(1,3,2)
    pl.plot(grid_x, phys_w)
    pl.title('w')
    pl.subplot(1,3,3)
    pl.plot(grid_x, phys_t)
    pl.title('T')
    pl.show()
