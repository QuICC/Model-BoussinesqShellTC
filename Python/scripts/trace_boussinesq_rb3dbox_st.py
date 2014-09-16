"""Script to run a marginal curve trace for the Boussinesq Rayleigh-Benard convection in a 3D box (streamfunction formulation)"""

import numpy as np

import geomhdiscc.model.boussinesq_rb3dbox_st as mod

# Create the model and activate linearization
model = mod.BoussinesqRB3DBoxST()
model.linearize = True
model.use_galerkin = False
fields = model.stability_fields()

# Set resolution, parameters, boundary conditions
res = [30, 30, 14]
#eq_params = {'prandtl':1, 'rayleigh':5555.993, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # Sol for NS/NS/NS, FF/FF/FT
#eq_params = {'prandtl':1, 'rayleigh':1886.5, 'scale1d':1.0/(12.0**0.5), 'scale2d':1.0/(12.0**0.5), 'scale3d':1.0} # Sol for NS/NS/NS, FF/FF/FT
eq_params = {'prandtl':1, 'rayleigh':1780.0, 'scale1d':1.0/5.0, 'scale2d':1.0/5.0, 'scale3d':1.0} # Sol for NS/NS/NS, FF/FF/FT
#eq_params = {'prandtl':1, 'rayleigh':3542.97, 'scale1d':1.0, 'scale2d':1.0, 'scale3d':1.0} # Sol for NS/NS/NS, FT/FT/FT
eigs = []
bc_str = 0 # 0: NS/NS/NS, 4: SF/SF/NS
bc_temp = 4 # 0: FT/FT/FT, 4: FF/FF/FT

bcs = {'bcType':model.SOLVER_HAS_BC, 'streamfunction':bc_str, 'temperature':bc_temp}

# Generate the operator A for the generalized EVP Ax = sigm B x
print("Constructing matrix A")
A = model.implicit_linear(res, eq_params, eigs, bcs, fields)

# Generate the operator B for the generalized EVP Ax = sigm B x
print("Constructing matrix B")
bcs['bcType'] = model.SOLVER_NO_TAU
B = model.time(res, eq_params, eigs, bcs, fields)

# Setup visualization and IO
show_spy = False
write_mtx = True
solve_evp = True
show_solution = (True and solve_evp)

if show_spy or show_solution:
    import matplotlib.pylab as pl
    #import mayavi.mbab as ml

if show_solution:
    import geomhdiscc.transform.cartesian as transf

# Show the "spy" of the two matrices
if show_spy:
    pl.spy(A, markersize=0.2)
    pl.show()
    pl.spy(B, markersize=0.2)
    pl.show()

# Export the two matrices to matrix market format
if write_mtx:
    import scipy.io as io
    print("Writing matrices to file")
    io.mmwrite("matrix_A.mtx", A)
    io.mmwrite("matrix_B.mtx", B)

# Solve EVP with sptarn
if solve_evp:
    print("Solve EVP")
    import geomhdiscc.linear_stability.solver as solver
    evp_vec, evp_lmb, iresult = solver.sptarn(A, B, -1e0, 10.0)
    print(evp_lmb)

if show_solution:
    viz_mode = -1
    xscale = eq_params['scale1d']
    yscale = eq_params['scale2d']
    zscale = eq_params['scale3d']

    print("\nVisualizing mode: " + str(evp_lmb[viz_mode]))
    # Get solution vectors
    sol_s = evp_vec[0:res[0]*res[1]*res[2],viz_mode]
    sol_t = evp_vec[res[0]*res[1]*res[2]:2*res[0]*res[1]*res[2],viz_mode]
    
    # Create spectrum plots
    pl.subplot(1,2,1)
    pl.semilogy(np.abs(sol_s))
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.semilogy(np.abs(sol_t))
    pl.title("T")
    pl.show()
    pl.close("all")
    
    # Create solution matrices
    mat_s = sol_s.reshape(res[0], res[1], res[2], order = 'F')
    mat_t = sol_t.reshape(res[0], res[1], res[2], order = 'F')

    # Compute physical space values
    grid_x = transf.grid(res[0])
    grid_y = transf.grid(res[1])
    grid_z = transf.grid(res[2])
    phys_s = transf.tophys3d(mat_s)
    phys_t = transf.tophys3d(mat_t)

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_z, grid_y, phys_s[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_z, grid_y, phys_t[res[0]//2,:,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_z, grid_x, phys_s[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_z, grid_x, phys_t[:,res[1]//2,:], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")

    # Show physical contour plot
    pl.subplot(1,2,1)
    pl.contourf(grid_y, grid_x, phys_s[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("Streamfunction")
    pl.subplot(1,2,2)
    pl.contourf(grid_y, grid_x, phys_t[:,:,res[2]//2], 50)
    pl.colorbar()
    pl.title("T")
    pl.show()
    pl.close("all")
