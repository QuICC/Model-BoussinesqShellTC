"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (velocity-continuity-temperature)"""

import numpy as np
import scipy.sparse as sp

import quicc.model.boussinesq_rotconv2dboxvct as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConv2DBoxVCT()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [60, 0, 60]
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':2340.687}
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.73}
eigs = [1]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityz':bc_vel, 'pressure':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
if False:
    import matplotlib.pylab as pl
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if True:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if True:
    import quicc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

    n_u = res[0]*(res[2]-1)
    n_w = (res[0]-1)*res[2]
    n_t = res[0]*res[2]
    n_p = (res[0]-1)*(res[2]-1)
    sol_u = evp_vec[0:n_u,-1]
    sol_w = evp_vec[n_u:n_u+n_w,-1]
    sol_t = evp_vec[n_u+n_w:n_u+n_w+n_t,-1]
    sol_p = evp_vec[n_u+n_w+n_t:n_u+n_w+n_t+n_p,-1]
    sol_c = mod.c2d.d1d0(res[0], res[2], {'x':{0:0}, 'z':{0:0, 'cr':1}}, sz = 0)*sol_u + mod.c2d.d0d1(res[0], res[2], {'x':{0:0, 'cr':1}, 'z':{0:0}}, sx = 0)*sol_w

    mat_u = (mod.c2d.qid(res[0], res[2], 0, 0, {'x':{0:0}, 'z':{0:0, 'cr':1}})*sol_u).reshape(res[0], res[2], order = 'F')
    mat_w = (mod.c2d.qid(res[0], res[2], 0, 0, {'x':{0:0, 'cr':1}, 'z':{0:0}})*sol_w).reshape(res[0], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[2], order = 'F')
    mat_p = (mod.c2d.qid(res[0], res[2], 0, 0, {'x':{0:0, 'cr':1}, 'z':{0:0, 'cr':1}})*sol_p).reshape(res[0], res[2], order = 'F')
    mat_c = sol_c.reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(2,3,1)
    pl.semilogy(np.abs(sol_u))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.semilogy(np.abs(sol_w))
    pl.title("w")
    pl.subplot(2,3,3)
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

    import matplotlib.pylab as pl
    pl.subplot(2,3,1)
    pl.imshow(np.log10(np.abs(mat_u)))
    pl.title("u")
    pl.subplot(2,3,2)
    pl.imshow(np.log10(np.abs(mat_w)))
    pl.title("w")
    pl.subplot(2,3,3)
    pl.imshow(np.log10(np.abs(mat_t)))
    pl.title("T")
    pl.subplot(2,3,5)
    pl.imshow(np.log10(np.abs(mat_c)))
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.imshow(np.log10(np.abs(mat_p)))
    pl.title("p")
    pl.show()
    pl.close("all")

    import quicc.transform.cartesian as transf
    phys_u = transf.tophys2d(mat_u)
    phys_w = transf.tophys2d(mat_w)
    phys_t = transf.tophys2d(mat_t)
    phys_p = transf.tophys2d(mat_p)
    phys_c = transf.tophys2d(mat_c)
    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(2,3,1)
    pl.contourf(grid_x, grid_z, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,3,2)
    pl.contourf(grid_x, grid_z, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,3,3)
    pl.contourf(grid_x, grid_z, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.subplot(2,3,5)
    pl.contourf(grid_x, grid_z, np.log10(np.abs(phys_c)), 50)
    pl.colorbar()
    pl.title("Continuity")
    pl.subplot(2,3,6)
    pl.contourf(grid_x, grid_z, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.show()
    pl.close("all")
