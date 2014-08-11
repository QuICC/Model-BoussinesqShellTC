"""Script to run a marginal curve trace for the Boussinesq rotating convection in a 2D box model (streamfunction-temperature)"""

import numpy as np

import geomhdiscc.model.boussinesq_rotconv2dboxst as mod

# Create the model and activate linearization
model = mod.BoussinesqRotConv2DBoxST()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]
eq_params = {'taylor':0, 'prandtl':1, 'rayleigh':5011.71461}
eigs = [1]
bc_str = 0 # 0: NS/NS + vel, 1: SF/SF + vel, 2: NS/NS + pressure, 3: SF/SF + pressure
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':bc_str, 'temperature':bc_temp}

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
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
    print(evp_lmb)

    sol_s = evp_vec[0:res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_t = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(1,2,1)
    pl.imshow(np.log10(np.abs(sol_s)))
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(np.log10(np.abs(sol_t)))
    pl.colorbar()
    pl.show()
    pl.close("all")

    import geomhdiscc.transform.cartesian as transf
    phys_s = transf.tophys2d(sol_s)
    phys_t = transf.tophys2d(sol_t)
    pl.subplot(1,2,1)
    pl.imshow(phys_s)
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.imshow(phys_t)
    pl.colorbar()
    pl.show()
    pl.close("all")

    grid_x = transf.grid(res[0])
    grid_z = transf.grid(res[2])

    pl.subplot(1,2,1)
    pl.contourf(grid_x, grid_z, phys_s, 50)
    pl.colorbar()
    pl.subplot(1,2,2)
    pl.contourf(grid_x, grid_z, phys_t, 50)
    pl.colorbar()
    pl.show()
    pl.close("all")
