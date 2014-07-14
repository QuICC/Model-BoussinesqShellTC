"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import numpy as np
import geomhdiscc.model.boussinesq_fplane3dqg as mod

# Create the model and activate linearization
model = mod.BoussinesqFPlane3DQG()
model.linearize = True
model.use_galerkin = True
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [50, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':6.6648, 'theta':35}

# Set wave number
phi = 0
kp = 1.2209
kx = kp*np.cos(phi*np.pi/180.0);
ky = (kp**2-kx**2)**0.5;
eigs = [kx, ky]

bcs = {'bcType':0, 'streamfunction':0, 'velocityz':0, 'temperature':0}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = 1
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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, 1)
    print(evp_lmb)
