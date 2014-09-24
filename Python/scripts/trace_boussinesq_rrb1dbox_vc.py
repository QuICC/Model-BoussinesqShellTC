"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a 1D box (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rrb1dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRRB1DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]

# SF, FT,
bc_vel = 1
bc_temp = 0
kx = 0
ky = 3.710
eq_params = {'prandtl':1, 'rayleigh':1676.12, 'taylor':1e3, 'scale1d':1.0}

eigs = [kx, ky]

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

# Show the "spy" of the two matrices
if False:
    import matplotlib.pylab as pl
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, np.inf)
    print(evp_lmb)

if show_solution:
    viz_mode = 0
    zscale = eq_params['scale1d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0],mode]
        sol_v = evp_vec[res[0]:2*res[0],mode]
        sol_w = evp_vec[2*res[0]:3*res[0],mode]

        # Extract continuity from velocity
        sol_c = 1j*(kx/2.0)*sol_u + 1j*(ky/2.0)*sol_v + mod.c1d.d1(res[0], mod.no_bc(), zscale)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0],viz_mode]
    sol_v = evp_vec[res[0]:2*res[0],viz_mode]
    sol_w = evp_vec[2*res[0]:3*res[0],viz_mode]
    sol_t = evp_vec[3*res[0]:4*res[0],viz_mode]
    sol_p = evp_vec[4*res[0]:5*res[0],viz_mode]
    # Extract continuity from velocity
    sol_c = 1j*(kx/2.0)*sol_u + 1j*(ky/2.0)*sol_v + mod.c1d.d1(res[0], mod.no_bc(), zscale)*sol_w

    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(abs(sol_u))
    pl.title('u')
    pl.subplot(2,3,2)
    pl.semilogy(abs(sol_v))
    pl.title('v')
    pl.subplot(2,3,3)
    pl.semilogy(abs(sol_w))
    pl.title('w')
    pl.subplot(2,3,4)
    pl.semilogy(abs(sol_t))
    pl.title('T')
    pl.subplot(2,3,5)
    pl.semilogy(abs(sol_p))
    pl.title('p')
    pl.subplot(2,3,6)
    pl.semilogy(abs(sol_c))
    pl.title('Continuity')
    pl.show()

    # Compute physical space values
    grid_x = transf.grid(res[0])
    phys_u = transf.tophys(sol_u.real)
    phys_v = transf.tophys(sol_v.real)
    phys_w = transf.tophys(sol_w.real)
    phys_t = transf.tophys(sol_t.real)
    phys_cr = transf.tophys(sol_c.real)
    phys_ci = transf.tophys(sol_c.imag)
    
    # Show physical plot
    pl.subplot(2,3,1)
    pl.plot(grid_x, phys_u)
    pl.title('u')
    pl.subplot(2,3,2)
    pl.plot(grid_x, phys_v)
    pl.title('v')
    pl.subplot(2,3,3)
    pl.plot(grid_x, phys_w)
    pl.title('w')
    pl.subplot(2,3,4)
    pl.plot(grid_x, phys_t)
    pl.title('T')
    pl.subplot(2,3,5)
    pl.plot(grid_x, phys_cr)
    pl.title('Continuity (real)')
    pl.subplot(2,3,6)
    pl.plot(grid_x, phys_ci)
    pl.title('Continuity (imag)')
    pl.show()