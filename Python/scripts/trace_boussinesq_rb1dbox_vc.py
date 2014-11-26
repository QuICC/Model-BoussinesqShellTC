"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a 1D box (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rb1dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRB1DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [64, 0, 0]

# SF, FT,
bc_vel = 1
bc_temp = 0
## kx = 0, ky = 2
#kx = 0
#ky = 2
#eq_params = {'prandtl':1, 'rayleigh':667.0098243, 'scale1d':2.0}
# kx = 0, ky = 1
#kx = 0
#ky = 1
#eq_params = {'prandtl':1, 'rayleigh':1284.225280, 'scale1d':2.0}
## kx = 0, ky = 5.35
#kx = 0
#ky = 5.35
#eq_params = {'prandtl':1, 'rayleigh':1992.541617, 'scale1d':2.0}
## kx = 2, ky = 0
#kx = 2
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':667.0098243, 'scale1d':2.0}
## kx = 1, ky = 0
#kx = 1
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':1284.225280, 'scale1d':2.0}
## kx = 5.35, ky = 0
#kx = 5.35
#ky = 0
#eq_params = {'prandtl':1, 'rayleigh':1992.541617, 'scale1d':2.0}
## Minimum n = 1: k_ = 2.221441469
#phi = 35
#kp = 2.221441469
#kx = kp*np.cos(phi*np.pi/180.0);
#ky = (kp**2-kx**2)**0.5;
#eq_params = {'prandtl':1, 'rayleigh':657.5113645, 'scale1d':2.0}
## Minimum n = 2: k_ = 4.442882938
#phi = 15
#kp = 4.442882938
#kx = kp*np.cos(phi*np.pi/180.0);
#ky = (kp**2-kx**2)**0.5;
#eq_params = {'prandtl':1, 'rayleigh':10520.18183, 'scale1d':2.0}
#kx = 1.31
#ky = 1.52
#eq_params = {'prandtl':1, 'rayleigh':833.12, 'scale1d':2.0}
# Minimum n = 2: k_ = 4.442882938
phi = 35
kp = np.sqrt(2)*np.pi
kx = kp*np.cos(phi*np.pi/180.0);
ky = (kp**2-kx**2)**0.5;
eq_params = {'prandtl':1, 'rayleigh':108*np.pi**4, 'scale1d':2.0}

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
show_solution = (False and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

# Show the "spy" of the two matrices
if show_spy:
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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -80.1e0, np.inf)
    print(evp_lmb)

    zscale = eq_params['scale1d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0],mode]
        sol_v = evp_vec[res[0]:2*res[0],mode]
        sol_w = evp_vec[2*res[0]:3*res[0],mode]

        # Extract continuity from velocity
        sol_c = 1j*kx*sol_u + 1j*ky*sol_v + mod.c1d.d1(res[0], mod.no_bc(), zscale)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

if show_solution:
    viz_mode = 0
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0],viz_mode]
    sol_v = evp_vec[res[0]:2*res[0],viz_mode]
    sol_w = evp_vec[2*res[0]:3*res[0],viz_mode]
    sol_t = evp_vec[3*res[0]:4*res[0],viz_mode]
    sol_p = evp_vec[4*res[0]:5*res[0],viz_mode]
    # Extract continuity from velocity
    sol_c = 1j*kx*sol_u + 1j*ky*sol_v + mod.c1d.d1(res[0], mod.no_bc(), zscale)*sol_w

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
