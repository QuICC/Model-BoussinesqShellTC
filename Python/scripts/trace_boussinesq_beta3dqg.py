"""Script to run a marginal curve trace for the Boussinesq Beta 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_beta3dqg as mod

# Create the model and activate linearization
model = mod.BoussinesqBeta3DQG()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [30, 0, 30]
chi = 1
eq_params = {'prandtl':1, 'rayleigh':700.0, 'gamma':1, 'chi':chi}
eigs = [2.221403788]

# No-slip/No-slip
#bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'temperature':0, 'velocityz':1} 
# Stress-free/Stress-free
#bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':1, 'temperature':0, 'velocityz':1}
# Stress-free/No-slip (simple Beta)
bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':2, 'temperature':0, 'velocityz':2, 'vorticityz':2}

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
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, 1)
    print(evp_lmb)

    sol_psi = evp_vec[0:res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_w = evp_vec[res[0]*res[2]:2*res[0]*res[2],-1].reshape(res[0], res[2], order = 'F')
    sol_t = evp_vec[2*res[0]*res[2]:,-1].reshape(res[0], res[2], order = 'F')

    import matplotlib.pylab as pl
    pl.subplot(1,3,1)
    pl.imshow(np.log10(np.abs(sol_psi)))
    pl.colorbar()
    pl.subplot(1,3,2)
    pl.imshow(np.log10(np.abs(sol_w)))
    pl.colorbar()
    pl.subplot(1,3,3)
    pl.imshow(np.log10(np.abs(sol_t)))
    pl.colorbar()
    pl.show()
    pl.close("all")
    import geomhdiscc.transform.cartesian as transf
    phys_psi = transf.tophys2d(sol_psi)
    phys_w = transf.tophys2d(sol_w)
    phys_t = transf.tophys2d(sol_t)
    pl.subplot(1,3,1)
    pl.imshow(sol_psi.real)
    pl.colorbar()
    pl.subplot(1,3,2)
    pl.imshow(sol_w.real)
    pl.colorbar()
    pl.subplot(1,3,3)
    pl.imshow(sol_t.real)
    pl.colorbar()
    pl.show()
    pl.close("all")

    xg = transf.grid(res[0])
    zg = transf.grid(res[2])
    pl.subplot(1,3,1)
    pl.plot(xg, sol_psi.real[:,10])
    pl.title('Psi')
    pl.subplot(1,3,2)
    pl.plot(xg, sol_w.real[:,10])
    pl.title('W')
    pl.subplot(1,3,3)
    pl.plot(xg, sol_t.real[:,10])
    pl.title('T')
    pl.show()
    pl.close("all")
    pl.subplot(1,3,1)
    pl.plot(zg, sol_psi.real[0,:])
    pl.title('Psi')
    pl.subplot(1,3,2)
    pl.plot(zg, sol_w.real[0,:])
    pl.title('W')
    pl.subplot(1,3,3)
    pl.plot(zg, sol_t.real[0,:])
    pl.title('T')
    pl.show()
    pl.close("all")
    pl.subplot(1,3,1)
    pl.plot(zg, sol_psi.real[res[0]//2,:])
    pl.title('Psi')
    pl.subplot(1,3,2)
    pl.plot(zg, sol_w.real[res[0]//2,:])
    pl.title('W')
    pl.subplot(1,3,3)
    pl.plot(zg, sol_t.real[res[0]//2,:])
    pl.title('T')
    pl.show()
    pl.close("all")
    pl.subplot(1,3,1)
    pl.plot(zg, sol_psi.real[-1,:])
    pl.title('Psi')
    pl.subplot(1,3,2)
    pl.plot(zg, sol_w.real[-1,:])
    pl.title('W')
    pl.subplot(1,3,3)
    pl.plot(zg, sol_t.real[-1,:])
    pl.title('T')
    pl.show()
    pl.close("all")
