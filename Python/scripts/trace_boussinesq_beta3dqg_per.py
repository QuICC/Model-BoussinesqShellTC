"""Script to run a marginal curve trace for the periodic Boussinesq Beta 3DQG model"""

import numpy as np
import geomhdiscc.model.boussinesq_beta3dqg_per as mod

# Create the model and activate linearization
model = mod.BoussinesqBeta3DQGPer()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [60, 0, 0]
chi = 60
Pr = 1
G = 1e-3
tanchi = np.tan(chi*np.pi/180.)
beta = tanchi/G

kxc = 0.0
kyc = (2.0**(-1.0/6.0)*(beta*Pr/(1.0 + Pr))**(1.0/3.0))
Rac = 16*((kxc**2 + kyc**2)**3/kyc**2 + beta**2*Pr**2/((kxc**2 + kyc**2)*(Pr + 1.0)**2))
fc = -beta*kyc/((kxc**2 + kyc**2)*(1.0 + Pr))
print("Ra_c = " + str(Rac))
print("k_c = " + str(kyc))
print("f_c = " + str(fc))

# Set wave number
kx = 0
ky = 14.4046937485972/2.0
eigs = [kx, ky]
Ra = 192407.5882

eq_params = {'prandtl':Pr, 'rayleigh':Ra, 'gamma':G, 'chi':chi, 'scale1d':1}
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0, 'vorticityz':0}

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
