"""Script to run a marginal curve trace for the Boussinesq rotating thermal convection in a spherical shell (Toroidal/Poloidal formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rtcshell as mod

# Create the model and activate linearization
model = mod.BoussinesqRTCShell()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
l = 0
m = 9
res = [20, m+20, 0]
eigs = [l, m]
eq_params = {'taylor':1e2, 'prandtl':1, 'rayleigh':4.761e6, 'ro':1, 'rratio':0.35}
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = True
write_mtx = True
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

if show_solution:
    viz_mode = 0
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    nL = res[1]-m+1
    sol_tor = evp_vec[0:res[0]*nL,viz_mode]
    sol_pol = evp_vec[res[0]*nL:2*res[0]*nL,viz_mode]
    sol_t = evp_vec[2*res[0]*nL:3*res[0]*nL,viz_mode]
    
    # Create spectrum plots
    pl.subplot(1,3,1)
    pl.semilogy(np.abs(sol_tor))
    pl.title("Toroidal")
    pl.subplot(1,3,2)
    pl.semilogy(np.abs(sol_pol))
    pl.title("Poloidal")
    pl.subplot(1,3,3)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.show()
    pl.close("all")
