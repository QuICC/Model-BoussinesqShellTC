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

set(GEOMHDISCC_FFTS "FFTW" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)

###################################################
#-------------- AVAILABLE FFT PLANS --------------#
###################################################

set(GEOMHDISCC_FFTPLANS "Fast" "Medium" "Slow" PARENT_SCOPE)

###################################################
#------- AVAILABLE LINEAR ALGEBRA LIBRARIES ------#
###################################################

set(GEOMHDISCC_LINALGS "LAPACK" "Eigen" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_LAPACK "auto" PARENT_SCOPE)

###################################################
#--- AVAILABLE SPARSE LINEAR ALGEBRA LIBRARIES ---#
###################################################

set(GEOMHDISCC_SPLINALGS "SuperLU" "UmfPack" "SparseLU" "KentLU" "Pardiso" "SPQR" "SparseQR" "BiCGSTAB" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SUPERLU "superlu" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SUPERLU "/usr/include/superlu" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_KENTLU "klu" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SPARSELU "/usr/local/lib/libmetis.a" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SPARSELU "/usr/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_PARDISO "/usr/local/Pardiso-4.1.2/lib" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#-------- AVAILABLE SPARSE EIGEN SOLVERS ---------#
###################################################

set(GEOMHDISCC_SPEIGSOLVERS "ARPACK" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_ARPACK "arpack" PARENT_SCOPE)

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

set(GEOMHDISCC_CC_ARCH_GCC "-march=native -O2" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_GCC "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_GCC ${GEOMHDISCC_CC_INC_GCC} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_GCC "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_GCC ${GEOMHDISCC_CC_LIB_GCC} PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_LIBRARIES "" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES "" PARENT_SCOPE)
