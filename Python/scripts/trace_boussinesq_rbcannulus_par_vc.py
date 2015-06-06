"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a cylindrical annulus model modified for parity (velocity-continuity formulation)"""

import numpy as np
import scipy.sparse as spsp
import scipy.sparse.linalg as spsplin

import geomhdiscc.model.boussinesq_rbcannulus_par_vc as mod

# Create the model and activate linearization
model = mod.BoussinesqRBCAnnulusVC()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
m = 2
res = [16, 0, 16]
#eq_params = {'prandtl':1, 'rayleigh':20412, 'ro':1, 'rratio':0.35, 'scale3d':2.0}
eq_params = {'prandtl':1, 'rayleigh':5e3, 'ro':1, 'rratio':0.35, 'scale3d':2.0}
eigs = [float(m)]
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 2 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = True
write_mtx = False
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.annulus as transf

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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e1, 1e2)
    print(evp_lmb)
    zscale = eq_params['scale3d']

    for mode in range(0,len(evp_lmb)):
        # Get solution vectors
        sol_ubar = evp_vec[0:res[0]*res[2],mode]
        sol_vbar = evp_vec[res[0]*res[2]:2*res[0]*res[2],mode]
        sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],mode]
        # Extract continuity from velocity 
        a, b = mod.annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])
        sol_c = mod.annulus.x1d1(res[0], res[2], a, b, mod.no_bc(), sr = 0, sz = 0)*sol_ubar + 1j*eigs[0]*sol_vbar + mod.annulus.x2e1(res[0], res[2], a, b, mod.no_bc(), zscale = zscale, sr = 0)*sol_w
        intg_c = mod.annulus.i1j1x1d1(res[0], res[2], a, b, mod.no_bc())*sol_ubar + mod.annulus.i1j1(res[0], res[2], a, b, mod.no_bc(),1j*eigs[0])*sol_vbar + mod.annulus.i1j1x2e1(res[0], res[2], a, b, mod.no_bc(), zscale = zscale)*sol_w
        print("Eigenvalue: " + str(evp_lmb[mode]) + ", Max continuity: " + str(np.max(np.abs(sol_c))) + ", Max integrated continuity: " + str(np.max(np.abs(intg_c))))

if show_solution:
    viz_mode = 0
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_ubar = evp_vec[0:res[0]*res[2],viz_mode]
    sol_vbar = evp_vec[res[0]*res[2]:2*res[0]*res[2],viz_mode]
    sol_w = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],viz_mode]
    sol_t = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],viz_mode]
    sol_pbar = evp_vec[4*res[0]*res[2]:5*res[0]*res[2],viz_mode]
    # Extract continuity from velocity 
    a, b = mod.annulus.rad.linear_r2x(eq_params['ro'], eq_params['rratio'])
    sol_c = mod.annulus.x1d1(res[0], res[2], a, b, mod.no_bc(), sr = 0, sz = 0)*sol_ubar + 1j*eigs[0]*sol_vbar + mod.annulus.x2e1(res[0], res[2], a, b, mod.no_bc(), zscale = zscale, sr = 0)*sol_w
    intg_c = mod.annulus.i1j1x1d1(res[0], res[2], a, b, mod.no_bc())*sol_ubar + mod.annulus.i1j1(res[0], res[2], a, b, mod.no_bc(),1j*eigs[0])*sol_vbar + mod.annulus.i1j1x2e1(res[0], res[2], a, b, mod.no_bc(), zscale = zscale)*sol_w
    
    # Create spectrum plots
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_ubar))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.semilogy(np.abs(sol_vbar))
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
    pl.semilogy(np.abs(sol_pbar))
    pl.title("p")
    pl.show()
    pl.close("all")
    
    # Create solution matrices
    mat_ubar = sol_ubar.reshape(res[0], res[2], order = 'F')
    mat_u = spsplin.spsolve(mod.annulus.rad.x1(res[0], a, b, {0:0}, zr = 0).tocsr(), mat_ubar.real)
    mat_vbar = sol_vbar.reshape(res[0], res[2], order = 'F')
    mat_v = spsplin.spsolve(mod.annulus.rad.x1(res[0], a, b, {0:0}, zr = 0).tocsr(), mat_vbar.real)
    mat_w = sol_w.reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_pbar = sol_pbar.reshape(res[0], res[2], order = 'F')
    mat_p = spsplin.spsolve(mod.annulus.rad.x2(res[0], a, b, {0:0}, zr = 0).tocsr(), mat_pbar.real)
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')

    # Visualize spectrum matrix
    pl.subplot(2,5,1)
    pl.imshow(np.log10(np.abs(mat_ubar)))
    pl.colorbar()
    pl.title("U bar")
    pl.subplot(2,5,2)
    pl.imshow(np.log10(np.abs(mat_vbar)))
    pl.colorbar()
    pl.title("V bar")
    pl.subplot(2,5,4)
    pl.imshow(np.log10(np.abs(mat_pbar)))
    pl.colorbar()
    pl.title("P bar")
    pl.subplot(2,5,5)
    pl.imshow(np.log10(np.abs(mat_c)))
    pl.colorbar()
    pl.title("Continuity")

    pl.subplot(2,5,6)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,5,7)
    pl.imshow(np.log10(np.abs(mat_v)))
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,5,8)
    pl.imshow(np.log10(np.abs(mat_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,5,9)
    pl.imshow(np.log10(np.abs(mat_p)))
    pl.colorbar()
    pl.title("p")
    pl.subplot(2,5,10)
    pl.imshow(np.log10(np.abs(mat_t)))
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    # Compute physical space values
    grid_r = transf.rgrid(res[0], a, b)
    grid_z = transf.zgrid(res[2])
    phys_ubar = transf.tophys2d(mat_ubar)
    phys_u = transf.tophys2d(mat_u)
    phys_vbar = transf.tophys2d(mat_vbar)
    phys_v = transf.tophys2d(mat_v)
    phys_w = transf.tophys2d(mat_w)
    phys_t = transf.tophys2d(mat_t)
    phys_pbar = transf.tophys2d(mat_pbar)
    phys_p = transf.tophys2d(mat_p)
    phys_c = transf.tophys2d(mat_c)

    # Show physical contour plot
    pl.subplot(2,5,1)
    pl.contourf(grid_z, grid_r, phys_ubar, 50)
    pl.colorbar()
    pl.title("U bar")
    pl.subplot(2,5,2)
    pl.contourf(grid_z, grid_r, phys_vbar, 50)
    pl.colorbar()
    pl.title("V bar")
    pl.subplot(2,5,4)
    pl.contourf(grid_z, grid_r, phys_pbar, 50)
    pl.colorbar()
    pl.title("P bar")
    pl.subplot(2,5,5)
    pl.contourf(grid_z, grid_r, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")

    pl.subplot(2,5,6)
    pl.contourf(grid_z, grid_r, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,5,7)
    pl.contourf(grid_z, grid_r, phys_v, 50)
    pl.colorbar()
    pl.title("v")
    pl.subplot(2,5,8)
    pl.contourf(grid_z, grid_r, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,5,9)
    pl.contourf(grid_z, grid_r, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.subplot(2,5,10)
    pl.contourf(grid_z, grid_r, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")