"""Script to run a marginal curve trace for the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation) without compled field coupling"""

import numpy as np

import geomhdiscc.model.boussinesq_tcshell_std as mod

# Create the model and activate linearization
model = mod.BoussinesqTCShellStd()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [20, 0, 0]
#l = 1; eq_params = {'prandtl':1, 'rayleigh':1.645e4, 'ro':1, 'rratio':0.2}
l = 1; eq_params = {'prandtl':1, 'rayleigh':3.0468e4, 'ro':1, 'rratio':0.3} #NS/NS
#l = 1; eq_params = {'prandtl':1, 'rayleigh':8.4791e3, 'ro':1, 'rratio':0.3} #SF/SF
#l = 1; eq_params = {'prandtl':1, 'rayleigh':4.183202e4, 'ro':1, 'rratio':0.5} #SF/SF
eigs = [l]
bc_vel = 0 # 0: NS/NS, 1: SF/SF, 2: SF/NS, 3: SF/NS
bc_temp = 0 # 0: FT/FT, 1: FF/FF, 2: FF/FT, 3: FT/FF
bcs = {'bcType':model.SOLVER_HAS_BC, 'velocity':bc_vel, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = False
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl

if show_solution:
    import geomhdiscc.transform.annulus as transf

# Show the "spy" of the two matrices
if show_spy:
    pl.spy(A, markersize=3, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()
    pl.spy(B, markersize=3, marker = '.', markeredgecolor = 'b')
    pl.tick_params(axis='x', labelsize=30)
    pl.tick_params(axis='y', labelsize=30)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e1, np.inf)
    print(evp_lmb)

if show_solution:
    viz_mode = 0
    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    if model.use_galerkin:
        sol_pol = evp_vec[0:res[0]-4,viz_mode]
        sol_t = evp_vec[res[0]-4:,viz_mode]
    else:
        sol_pol = evp_vec[0:res[0],viz_mode]
        sol_t = evp_vec[res[0]:2*res[0],viz_mode]
    
    # Create spectrum plots
    pl.subplot(1,2,1)
    pl.semilogy(np.abs(sol_pol))
    pl.title("Poloidal")
    pl.subplot(1,2,2)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.show()
    pl.close("all")
