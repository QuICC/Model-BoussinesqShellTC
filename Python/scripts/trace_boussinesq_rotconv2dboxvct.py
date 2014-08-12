"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (velocity-continuity-temperature)"""

import numpy as np

import geomhdiscc.model.boussinesq_rotconv2dboxvct as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConv2DBoxVCT()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [10, 0, 10]
#eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':2340.687}
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011}
eigs = [1]
bc_vel = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'velocityx':bc_vel, 'velocityz':bc_vel, 'pressure':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)
A = A.tolil()
#A[3*res[0]*res[2],:] = 0
#A[3*res[0]*res[2],3*res[0]*res[2]] = 1
#A[3*res[0]*res[2],3*res[0]*res[2]] = 1
#A[3*res[0]*res[2]+1,3*res[0]*res[2]+res[0]-1] = 1
#A[3*res[0]*res[2]+res[0],-res[0]] = 1
#A[3*res[0]*res[2]+res[0]+1,-1] = 1
A = A.tocsr()

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
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

    sol_u = evp_vec[0:res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_w = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_t = evp_vec[2*res[0]*res[2]:3*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_p = evp_vec[3*res[0]*res[2]:4*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(2,2,1)
    pl.imshow(np.log10(np.abs(sol_u)))
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,2,2)
    pl.imshow(np.log10(np.abs(sol_w)))
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,2,3)
    pl.imshow(np.log10(np.abs(sol_p)))
    pl.colorbar()
    pl.title("p")
    pl.subplot(2,2,4)
    pl.imshow(np.log10(np.abs(sol_t)))
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    import geomhdiscc.transform.cartesian as transf
    phys_u = transf.tophys2d(sol_u)
    phys_w = transf.tophys2d(sol_w)
    phys_t = transf.tophys2d(sol_t)
    phys_p = transf.tophys2d(sol_p)
    pl.subplot(2,2,1)
    pl.imshow(phys_u)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,2,2)
    pl.imshow(phys_w)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,2,3)
    pl.imshow(phys_p)
    pl.colorbar()
    pl.title("p")
    pl.subplot(2,2,4)
    pl.imshow(phys_t)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(2,2,1)
    pl.contourf(grid_x, grid_z, phys_u, 50)
    pl.colorbar()
    pl.title("u")
    pl.subplot(2,2,2)
    pl.contourf(grid_x, grid_z, phys_w, 50)
    pl.colorbar()
    pl.title("w")
    pl.subplot(2,2,3)
    pl.contourf(grid_x, grid_z, phys_p, 50)
    pl.colorbar()
    pl.title("p")
    pl.subplot(2,2,4)
    pl.contourf(grid_x, grid_z, phys_t, 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")
