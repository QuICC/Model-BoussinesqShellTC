"""Script to run a marginal curve trace for the Boussinesq Beta 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_beta3dqg2 as mod

# Create the model and activate linearization
model = mod.BoussinesqBeta3DQG()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 20]
chi = 1
eq_params = {'prandtl':1, 'rayleigh':1711.5, 'gamma':1, 'chi':chi}
eigs = [3.11627]
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'phi':0, 'temperature':0, 'velocityz':0, 'lambda':0, 'chi':0}

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
if False:
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, 1)
    print(evp_lmb)
