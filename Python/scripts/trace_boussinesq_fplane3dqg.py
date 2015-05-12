"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_fplane3dqg as mod

# Create the model and activate linearization
model = mod.BoussinesqFPlane3DQG()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'scale1d':2.0}
res = [32, 0, 0]
eq_params = {'prandtl':7, 'rayleigh':8.6957, 'scale1d':2.0}

# Set wave number
phi = 0
kp = 1.3048
kx = kp*np.cos(phi*np.pi/180.0);
ky = (kp**2-kx**2)**0.5;
eigs = [kx, ky]

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = False
solve_evp = True
show_solution = (False and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

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

# Solve EVP with sptarn
if solve_evp:
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
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
