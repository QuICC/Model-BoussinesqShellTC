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

set(GEOMHDISCC_SMARTPTRS "TR1" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(GEOMHDISCC_THREADSMODELS "None" PARENT_SCOPE)

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

set(GEOMHDISCC_LINALGS "Eigen" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_LAPACK "auto" PARENT_SCOPE)

###################################################
#--- AVAILABLE SPARSE LINEAR ALGEBRA LIBRARIES ---#
###################################################

set(GEOMHDISCC_SPLINALGS "UmfPack" "SparseLU" "SPQR" "SparseQR" "BiCGSTAB" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_UMFPACK "/usr/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "suitesparseconfig" "lapack" "blas" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_UMFPACK "/usr/local/lib" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SPARSELU "/usr/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(GEOMHDISCC_MPIIMPLS "OpenMPI" PARENT_SCOPE)

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
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_PYTHONS "python34" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON34 "/usr/lib/python3.4/config-3.4m-x86_64-linux-gnu/libpython3.4.so" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON34 "/usr/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_LIBRARIES "" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES "" PARENT_SCOPE)
