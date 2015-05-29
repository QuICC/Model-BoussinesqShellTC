"""Script to run a marginal curve trace for the Boussinesq F-plane 3DQG model"""

import numpy as np

import geomhdiscc.model.boussinesq_fplane3dqg as mod
import geomhdiscc.linear_stability.MarginalCurve as MarginalCurve

# Create the model and activate linearization
model = mod.BoussinesqFPlane3DQG()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [32, 0, 0]
eq_params = {'prandtl':1, 'rayleigh':8.6957, 'scale1d':2.0}
res = [32, 0, 0]
eq_params = {'prandtl':7, 'rayleigh':8.6957, 'scale1d':2.0}

# Set wave number
def wave(kp):
    phi = 0
    kx = kp*np.cos(phi*np.pi/180.0)
    ky = (kp**2-kx**2)**0.5
    return [kx, ky]

eigs = [1, 0]

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':0, 'velocityz':0, 'temperature':0}

gevp = MarginalCurve.GEVP(model, res, eq_params, eigs, bcs, fields, wave = wave)

curve = MarginalCurve.MarginalCurve(gevp)
ks = np.arange(1, 2, 0.05)

(data_k, data_Ra, data_freq) = curve(ks)
kc, Rac, fc = curve.minimum(data_k, data_Ra)
import matplotlib.pylab as pl
pl.plot(data_k, data_Ra, 'b', data_k, data_freq, 'g', kc, Rac, 'r+', markersize=12)
pl.show()

curve.point.gevp.setEigs(kc)
evp_vec = curve.point.eigenvector()
evp_lmb = curve.point.eigenvalue()
evp_lmb = curve.point.gevp.solve(2*Rac, 10)
print(curve.point.gevp.evp_lmb)

# Setup visualization and IO
show_spy = False
write_mtx = False
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.cartesian as transf

# Show the "spy" of the two matrices
if show_spy:
    A, B = curve.gevp(1.0)
    pl.spy(A, markersize=5, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()
    pl.spy(B, markersize=5, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
#if solve_evp:
#    import geomhdiscc.linear_stability.solver as solver
#    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1, np.inf)
#    print(evp_lmb)

if show_solution:
#    viz_mode = 0
#    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
#    sol_s = evp_vec[0:res[0],viz_mode]
#    sol_w = evp_vec[res[0]:2*res[0],viz_mode]
#    sol_t = evp_vec[2*res[0]:3*res[0],viz_mode]
    print("\nVisualizing mode: " + str(evp_lmb))
    sol_s = evp_vec[0:res[0]]
    sol_w = evp_vec[res[0]:2*res[0]]
    sol_t = evp_vec[2*res[0]:3*res[0]]
    # Create spectrum plots
    pl.subplot(1,3,1)
    pl.semilogy(abs(sol_s))
    pl.title('s')
    pl.subplot(1,3,2)
    pl.semilogy(abs(sol_w))
    pl.title('w')
    pl.subplot(1,3,3)
    pl.semilogy(abs(sol_t))
    pl.title('T')
    pl.show()

    # Compute physical space values
    grid_x = transf.grid(res[0])
    phys_s = transf.tophys(sol_s.real)
    phys_w = transf.tophys(sol_w.real)
    phys_t = transf.tophys(sol_t.real)

    # Show physical plot
    pl.subplot(1,3,1)
    pl.plot(grid_x, phys_s)
    pl.title('s')
    pl.subplot(1,3,2)
    pl.plot(grid_x, phys_w)
    pl.title('w')
    pl.subplot(1,3,3)
    pl.plot(grid_x, phys_t)
    pl.title('T')
    pl.show()
