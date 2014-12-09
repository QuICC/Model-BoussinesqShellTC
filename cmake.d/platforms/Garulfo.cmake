###################################################
#----------- Desktop CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(GEOMHDISCC_COMPILERS "GCC" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(GEOMHDISCC_SMARTPTRS "Boost" "TR1" "cxx0x" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(GEOMHDISCC_FFTS "FFTW" "cuFFT" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_CUFFT "cudart" "cufft" "cufftw" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_CUFFT "/opt/cuda/include" "/opt/cuda/sdk/common/inc/" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_CUFFT "/opt/cuda/lib64" PARENT_SCOPE)

###################################################
#-------------- AVAILABLE FFT PLANS --------------#
###################################################

set(GEOMHDISCC_FFTPLANS "Fast" "Medium" "Slow" PARENT_SCOPE)

###################################################
#------- AVAILABLE LINEAR ALGEBRA LIBRARIES ------#
###################################################

set(GEOMHDISCC_LINALGS "Eigen" PARENT_SCOPE)

###################################################
#--- AVAILABLE SPARSE LINEAR ALGEBRA LIBRARIES ---#
###################################################

set(GEOMHDISCC_SPLINALGS "SuperLU" "UmfPack" "SparseLU" "Mumps" "Pardiso" "SPQR" "SparseQR" "BiCGSTAB" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SUPERLU "superlu" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SUPERLU "/usr/include/superlu" PARENT_SCOPE)
#set(GEOMHDISCC_LIBRARIES_SPARSELU "/usr/local/lib/libmetis.a" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SPARSELU "/usr/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MUMPS "dmumps" "zmumps" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_PARDISO "pardiso412-GNU450-X86-64" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_PARDISO "/usr/local/share/Pardiso/lib" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#--------- AVAILABLE TRANSFORM OPERATORS ---------#
###################################################

set(GEOMHDISCC_TRANSOPS "Forward" "Recurrence" "Backward" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(GEOMHDISCC_MPLIBS "mpfr" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(GEOMHDISCC_CC_SERIAL_GCC "g++" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_GCC "mpic++" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_GCC "-march=native -O2 -fopenmp" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_GCC "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_GCC ${GEOMHDISCC_CC_INC_GCC} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_GCC "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_GCC ${GEOMHDISCC_CC_LIB_GCC} PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

# Find python headers and library
find_package(PythonLibs REQUIRED)

set(GEOMHDISCC_LIBRARIES ${PYTHON_LIBRARIES} PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES ${PYTHON_INCLUDE_DIRS} PARENT_SCOPE)
