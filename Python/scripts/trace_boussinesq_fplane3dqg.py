"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import geomhdiscc.model.boussinesq_fplane3dqg

# Create the model and activate linearization
model = geomhdiscc.model.boussinesq_fplane3dqg.BoussinesqFPlane3DQG()
model.linearize = True
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [50, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':4.761e6, 'theta':1}
eigs = [2, 3]
bcs = {'bcType':0, 'streamfunction':0, 'velocityz':0, 'temperature':0}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = 2
B = model.time(res, eq_params, eigs, bcs, fields)

# Show the "spy" of the two matrices
#import matplotlib.pylab as pl
#pl.spy(A, markersize=0.2)
#pl.show()
#pl.spy(B, markersize=0.2)
#pl.show()

# Export the two matrices to matrix market format
import scipy.io as io
io.mmwrite("matrix_A.mtx", A)
io.mmwrite("matrix_B.mtx", B)
