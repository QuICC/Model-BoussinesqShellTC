"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a 3D box (velocity-continuity formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rb3dbox_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRB3DBoxVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [14, 14, 14]
#eq_params = {'prandtl':1, 'rayleigh':2340.687, 'zxratio':1.0, 'yxratio':1.0}
eq_params = {'prandtl':1, 'rayleigh':13011.73, 'zxratio':1.0, 'yxratio':1.0}
eigs = []
bc_vel = 0 # 0: NS/NS/NS, 1: SF/SF/SF, 2: SF/NS/NS, 3: SF/NS/NS
bc_temp = 0 # 0: FT/FT/FT, 1: FF/FF/FF, 2: FF/FT/FT, 3: FT/FF/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityy':bc_vel, 'velocityz':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
print("Constructing matrix A")
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
print("Constructing matrix B")
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl
    #import mayavi.mbab as ml

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
    print("Solve EVP")
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

if show_solution:
    viz_mode = -1

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_u = evp_vec[0:res[0]*res[1]*res[2],mode]
        sol_v = evp_vec[res[0]*res[1]*res[2]:2*res[0]*res[1]*res[2],mode]
        sol_w = evp_vec[2*res[0]*res[1]*res[2]:3*res[0]*res[1]*res[2],mode]
        # Extract continuity from velocity 
        sol_c = mod.c3d.d1(res[0], res[1], res[2], mod.no_bc(), sy = 0, sz = 0)*sol_u + mod.c3d.e1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sz = 0)*sol_v + mod.c3d.f1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sy = 0)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))))

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_u = evp_vec[0:res[0]*res[1]*res[2],viz_mode]
    sol_v = evp_vec[res[0]*res[1]*res[2]:2*res[0]*res[1]*res[2],viz_mode]
    sol_w = evp_vec[2*res[0]*res[1]*res[2]:3*res[0]*res[1]*res[2],viz_mode]
    sol_t = evp_vec[3*res[0]*res[1]*res[2]:4*res[0]*res[1]*res[2],viz_mode]
    sol_p = evp_vec[4*res[0]*res[1]*res[2]:5*res[0]*res[1]*res[2],viz_mode]
    # Extract continuity from velocity 
    sol_c = mod.c3d.d1(res[0], res[1], res[2], mod.no_bc(), sy = 0, sz = 0)*sol_u + mod.c3d.e1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sz = 0)*sol_v + mod.c3d.f1(res[0], res[1], res[2], mod.no_bc(), sx = 0, sy = 0)*sol_w
    
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
    mat_u = sol_u.reshape(res[0], res[1], res[2], order = 'F')
    mat_v = sol_v.reshape(res[0], res[1], res[2], order = 'F')
    mat_w = sol_w.reshape(res[0], res[1], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[1], res[2], order = 'F')
    mat_p = sol_p.reshape(res[0], res[1], res[2], order = 'F')
    mat_c = sol_c.reshape(res[0], res[1], res[2], order = 'F')

    # Compute physical space values
    grid_x = transf.grid(res[0])
    grid_y = transf.grid(res[1])
    grid_z = transf.grid(res[2])
    phys_u = transf.tophys3d(mat_u)
    phys_v = transf.tophys3d(mat_v)
    phys_w = transf.tophys3d(mat_w)
    phys_t = transf.tophys3d(mat_t)
    phys_p = transf.tophys3d(mat_p)
    phys_c = transf.tophys3d(mat_c)

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_v[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_w[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_x, grid_z, phys_t[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c[res[0]//2,:,:])), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_v[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_w[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_x, grid_z, phys_t[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c[:,res[1]//2,:])), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_v[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_w[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,4)
    pl.contourf(grid_x, grid_z, phys_t[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c[:,:,res[2]//2])), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
