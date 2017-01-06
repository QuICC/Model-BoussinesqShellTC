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

set(QUICC_SMARTPTRS "Boost" "TR1" "cxx0x" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "SuperLU" "UmfPack" "SparseLU" "MUMPS" "SPQR" "SparseQR" "BiCGSTAB" PARENT_SCOPE)
set(QUICC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(QUICC_LIBRARIES_SUPERLU "superlu" PARENT_SCOPE)
set(QUICC_INCLUDES_SUPERLU "/usr/include/superlu" PARENT_SCOPE)
#set(QUICC_LIBRARIES_SPARSELU "/usr/local/lib/libmetis.a" PARENT_SCOPE)
set(QUICC_INCLUDES_SPARSELU "/usr/local/include" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" PARENT_SCOPE)
set(QUICC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "SuperLU" "UmfPack" "SparseLU" "MUMPS" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "SuperLU" "UmfPack" "MUMPS" PARENT_SCOPE)

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

set(QUICC_MPLIBS "mpfr" PARENT_SCOPE)
set(QUICC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(QUICC_CC_SERIAL_GCC "g++" PARENT_SCOPE)

set(QUICC_CC_MPI_GCC "mpic++" PARENT_SCOPE)

set(QUICC_CC_ARCH_GCC "-march=native -O2" PARENT_SCOPE)

set(QUICC_CC_INC_GCC "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_GCC ${QUICC_CC_INC_GCC} PARENT_SCOPE)

set(QUICC_CC_LIB_GCC "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_GCC ${QUICC_CC_LIB_GCC} PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" "python34" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON27 "/usr/lib64/libpython2.7.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27 "/usr/include/python2.7" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON34 "/usr/lib64/libpython3.4.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON34 "/usr/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
