# Quasi-Inverse Convection Code

### QuICC is a simulation framework for solving the Navier-Stokes equations in spherical, cylindrical and Cartesian geometry.

## Features:

   - Fully spectral using Fourier, Chebyshev, Worland and Associated Legendre basis
   - Multiple geometries: spherical shell, full sphere, cylinder, plane layer, cartesian box
   - Sparse formulation of the linear operators
   - High-level Python interface for the linear operators
   - Nonlinear simulations
   - Linear stability calculations using PETSc/SLEPc
   - 1D, 2D and 2.5D parallelisation

## Quick Start

0. Upload your SSH key to your gitlab profile
1. Clone repository
   ```
   git clone git@gitlab.base13.ch:QuICC/QuICC.git
   ```
2. Clone submodules
   ```
   cd /path/to/QuICC
   git submodule init
   git submodule update
   ```
3. Create build directory
   ```
   mkdir /path/to/Builds
   cd /path/to/Builds
   ```
4. Configure build with CMake
   ```
   ccmake /path/to/QuICC
   ```
5. Compile the model executables
   ```
   make GreatSimulationConfig
   make GreatSimulationState
   make GreatSimulationModel
   make GreatSimulationVisu
   ```
6. Create configuration XML
   ```
   ./GreatSimulationConfig
   mv parameters_GEOMETRY.cfg parameters.cfg
   edit parameters.cfg
   ```
7. Create initial state
   ```
   ./GreatSimulationState
   mv state0000.hdf5 state_initial.hdf5
   ```
8. Run simulation
   ```
   ./GreatSimulationModel
   ```
9. Create physical space data for visualization
   ```
   ln -s state0042.hdf5 state4Visu.hdf5
   ./GreatSimulationVisu
   ```
10. visualize *visState0000.hdf5*