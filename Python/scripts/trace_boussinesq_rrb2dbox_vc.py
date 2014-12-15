"""Script to run a marginal curve trace for the Boussinesq rotating Rayleigh-Benard convection in a 2D box (1 periodic direction) (velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rrb2dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRRB2DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [14, 0, 14]

## SF/SF, FF/FT, k = 0
bc_vel = 2
bc_temp = 1
eigs = [5]
## SF/SF, FF/FT, Aspect ratio 1:1
eq_params = {'prandtl':1, 'rayleigh':1e5, 'taylor':1e3, 'scale1d':2.0, 'scale3d':2.0} # m = 1, n = 1, aspect ration 1:1

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}


# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = True
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

# Show the "spy" of the two matrices
if show_spy:
    import matplotlib.pylab as pl
    pl.spy(A, markersize=3, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()
    pl.spy(B, markersize=3, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    print("Writing MTX files")
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, np.inf)
    print("Found " + str(len(evp_lmb)) + " eigenvalues\n")

if show_solution:
    viz_mode = 0
    k = eigs[0]
    xscale = eq_params['scale1d']
    zscale = eq_params['scale3d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0]*res[2],mode]
        sol_v = evp_vec[res[0]*res[2]:2*res[0]*res[2],mode]
        sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],mode]

        # Extract continuity from velocity 
        sol_c = mod.c2d.d1(res[0], res[2], mod.no_bc(), xscale = xscale, sz = 0)*sol_u + 1j*k*sol_v + mod.c2d.e1(res[0], res[2], mod.no_bc(), zscale = zscale, sx = 0)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0]*res[2],viz_mode]
    sol_v = evp_vec[res[0]*res[2]:2*res[0]*res[2],viz_mode]
    sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],viz_mode]
    sol_t = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],viz_mode]
    sol_p = evp_vec[4*res[0]*res[2]:5*res[0]*res[2],viz_mode]

    # Extract continuity from velocity 
    sol_c = mod.c2d.d1(res[0], res[2], mod.no_bc(), xscale = xscale, sz = 0)*sol_u + 1j*k*sol_v + mod.c2d.e1(res[0], res[2], mod.no_bc(), zscale = zscale, sx = 0)*sol_w
    
    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_u))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.semilogy(np.abs(sol_v))
    pl.title("v")
    pl.subplot(2,3,3)
    pl.semilogy(np.abs(sol_w))
    pl.title("w")
    pl.subplot(2,3,4)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.subplot(2,3,5)
    pl.semilogy(np.abs(sol_c))
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.semilogy(np.abs(sol_p))
    pl.title("p")
    pl.show()
    pl.close("all")
    
    # Create solution matrices
    mat_u = sol_u.reshape(res[0], res[2], order = 'F')
    mat_v = sol_v.reshape(res[0], res[2], order = 'F')
    mat_w = sol_w.reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_p = sol_p.reshape(res[0], res[2], order = 'F')
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')

    # Visualize spectrum matrix
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(mat_v)))
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(mat_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.imshow(np.log10(np.abs(mat_t)))
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(mat_c)))
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(mat_p)))
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    # Compute physical space values
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])
    phys_u = transf.tophys2d(mat_u)
    phys_v = transf.tophys2d(mat_v)
    phys_w = transf.tophys2d(mat_w)
    phys_t = transf.tophys2d(mat_t)
    phys_p = transf.tophys2d(mat_p)
    phys_c = transf.tophys2d(mat_c)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_z, grid_x, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_z, grid_x, phys_v, 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.contourf(grid_z, grid_x, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_z, grid_x, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_z, grid_x, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_z, grid_x, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
