"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a 1D box (2 periodic directions) (Velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rrb1dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRRB1DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [128, 0, 0]

# FT
bc_temp = 0

# SF, FT,
#bc_vel = 1
kx = 0
ky = 3.710
eq_params = {'prandtl':1, 'rayleigh':1676.12, 'taylor':1e3, 'heating':0, 'scale1d':2.0}
#ky = 129
#eq_params = {'prandtl':1, 'rayleigh':8.7050552e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 131
#eq_params = {'prandtl':1, 'rayleigh':8.7012659e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 130
#eq_params = {'prandtl':1, 'rayleigh':8.7011095e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 132
#eq_params = {'prandtl':1, 'rayleigh':8.7054933e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 65
#eq_params = {'prandtl':1, 'rayleigh':4.13360252e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60
#eq_params = {'prandtl':1, 'rayleigh':4.04824521e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 61
#eq_params = {'prandtl':1, 'rayleigh':4.04803724e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 62
#eq_params = {'prandtl':1, 'rayleigh':4.05657944e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60.5
#eq_params = {'prandtl':1, 'rayleigh':4.04703887e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60.52
#eq_params = {'prandtl':1, 'rayleigh':4.04703659e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60.523
#eq_params = {'prandtl':1, 'rayleigh':4.04703656e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60.5229
#eq_params = {'prandtl':1, 'rayleigh':4.04703656e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
## NS, FT,
bc_vel = 0
#ky = 54
#eq_params = {'prandtl':1, 'rayleigh':3.455838e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
ky = 55
eq_params = {'prandtl':1, 'rayleigh':3.450289e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
ky = 55.4
eq_params = {'prandtl':1, 'rayleigh':3.44979010e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
ky = 55.402
eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
ky = 55.4023
eq_params = {'prandtl':1, 'rayleigh':3.44979009e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 56
#eq_params = {'prandtl':1, 'rayleigh':3.450893e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 60
#eq_params = {'prandtl':1, 'rayleigh':3.515608e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 70
#eq_params = {'prandtl':1, 'rayleigh':4.134931e7, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 100
#eq_params = {'prandtl':1, 'rayleigh':1.094775e8, 'taylor':1e10, 'heating':0, 'scale1d':2.0}
#ky = 120
#eq_params = {'prandtl':1, 'rayleigh':7.7971364e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 121
#eq_params = {'prandtl':1, 'rayleigh':7.7892799e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 122
#eq_params = {'prandtl':1, 'rayleigh':7.7846430e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 123
#eq_params = {'prandtl':1, 'rayleigh':7.7832219e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 124
#eq_params = {'prandtl':1, 'rayleigh':7.7850143e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 128
#eq_params = {'prandtl':1, 'rayleigh':7.8420000e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 130
#eq_params = {'prandtl':1, 'rayleigh':7.8632272e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}
#ky = 131
#eq_params = {'prandtl':1, 'rayleigh':7.8875224e8, 'taylor':1e12, 'heating':0, 'scale1d':2.0}

eigs = [kx, ky]

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = False
solve_evp = True
show_solution = (True and solve_evp)

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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e3, np.inf)
    print(evp_lmb)

if show_solution:
    viz_mode = -1
    zscale = eq_params['scale1d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0],mode]
        sol_v = evp_vec[res[0]:2*res[0],mode]
        sol_w = evp_vec[2*res[0]:3*res[0],mode]

        # Extract continuity from velocity
        sol_c = 1j*kx*sol_u + 1j*ky*sol_v + mod.c1d.d1(res[0], mod.no_bc(), zscale)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

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
    phys_ur = transf.tophys(sol_u.real)
    phys_ui = transf.tophys(sol_u.imag)
    phys_vr = transf.tophys(sol_v.real)
    phys_vi = transf.tophys(sol_v.imag)
    phys_wr = transf.tophys(sol_w.real)
    phys_wi = transf.tophys(sol_w.imag)
    phys_tr = transf.tophys(sol_t.real)
    phys_ti = transf.tophys(sol_t.imag)
    phys_pr = transf.tophys(sol_p.real)
    phys_pi = transf.tophys(sol_p.imag)
    phys_cr = transf.tophys(sol_c.real)
    phys_ci = transf.tophys(sol_c.imag)
    
    # Show physical plot
    pl.subplot(2,3,1)
    pl.plot(grid_x, phys_ur)
    pl.plot(grid_x, phys_ui)
    pl.title('u')
    pl.subplot(2,3,2)
    pl.plot(grid_x, phys_vr)
    pl.plot(grid_x, phys_vi)
    pl.title('v')
    pl.subplot(2,3,3)
    pl.plot(grid_x, phys_wr)
    pl.plot(grid_x, phys_wi)
    pl.title('w')
    pl.subplot(2,3,4)
    pl.plot(grid_x, phys_tr)
    pl.plot(grid_x, phys_ti)
    pl.title('T')
    pl.subplot(2,3,5)
    pl.plot(grid_x, phys_pr)
    pl.plot(grid_x, phys_pi)
    pl.title('Pressure (real)')
    pl.subplot(2,3,6)
    pl.plot(grid_x, phys_ci)
    pl.plot(grid_x, phys_cr)
    pl.title('Continuity (imag)')
    pl.show()
