###################################################
#----------- Desktop CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "GCC" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(QUICC_SMARTPTRS "std" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(QUICC_THREADSMODELS "None" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(QUICC_FFTS "FFTW" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)

###################################################
#-------------- AVAILABLE FFT PLANS --------------#
###################################################

set(QUICC_FFTPLANS "Fast" "Medium" "Slow" PARENT_SCOPE)

###################################################
#------- AVAILABLE LINEAR ALGEBRA LIBRARIES ------#
###################################################

set(QUICC_LINALGS "Eigen" PARENT_SCOPE)

###################################################
#--- AVAILABLE SPARSE LINEAR ALGEBRA LIBRARIES ---#
###################################################

set(QUICC_SPLINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)
# MUMPS
#set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "openblas" "goto2" "mpi" "mpi_f90" "mpi_f77" "gfortran" "cmumps" "smumps" "pord" "scalapack" "metis" "parmetis" PARENT_SCOPE)
#set(QUICC_INCLUDES_MUMPS "$ENV{MUMPS_ROOT}/include" PARENT_SCOPE)
#set(QUICC_LIBDIR_MUMPS "$ENV{MUMPS_ROOT}/lib" "$ENV{MODULE_open_mpi_ROOT_DIR}/lib" PARENT_SCOPE)
# UmfPack
#set(QUICC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
#set(QUICC_INCLUDES_UMFPACK "$ENV{SUITE_SPARSE_ROOT}/include" PARENT_SCOPE)
#set(QUICC_LIBDIR_UMFPACK "$ENV{SUITE_SPARSE_ROOT}/lib" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "UmfPack" "SparseLU" "MUMPS" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(QUICC_MPIIMPLS "OpenMPI" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(QUICC_MPBACKENDS "boost" "gmp" "mpfr" "quad" PARENT_SCOPE)
set(QUICC_INCLUDES_BOOST "/usr/include" PARENT_SCOPE)
set(QUICC_INCLUDES_GMP "/usr/include" PARENT_SCOPE)
set(QUICC_LIBRARIES_GMP "gmp" PARENT_SCOPE)
set(QUICC_INCLUDES_MPFR "/usr/include" PARENT_SCOPE)
set(QUICC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)
set(QUICC_LIBRARIES_QUAD "quadmath" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(QUICC_CC_SERIAL_GCC "g++" PARENT_SCOPE)

set(QUICC_CC_MPI_GCC "mpic++" PARENT_SCOPE)

set(QUICC_CC_ARCH_GCC "-march=native -O2 -std=c++11" PARENT_SCOPE)

set(QUICC_CC_OPENMP_GCC "-fopenmp" PARENT_SCOPE)

set(QUICC_CC_INC_GCC "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_GCC "" PARENT_SCOPE)

set(QUICC_CC_LIB_GCC "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_GCC "" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python36" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON36 "$ENV{PYTHON_ROOT}/lib/libpython3.6m.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON36 "$ENV{PYTHON_ROOT}/include/python3.6m" "$ENV{PYTHON_ROOT}/lib/python3.6/site-packages/numpy/core/include" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
