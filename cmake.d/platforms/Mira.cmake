###################################################
#----------- Janus CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "Clang" PARENT_SCOPE)

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
set(QUICC_INCLUDES_FFTW_CLANG "/soft/libraries/alcf/current/gcc/FFTW3/include" PARENT_SCOPE)
set(QUICC_LIBDIR_FFTW_CLANG "/soft/libraries/alcf/current/gcc/FFTW3/lib" PARENT_SCOPE)

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
set(QUICC_LIBRARIES_UMFPACK "umfpack" "m" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" "lapack" "blas" "gfortran" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_CLANG "/soft/libraries/petsc/3.5.2.1/xl-opt/include/" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_CLANG "/soft/libraries/alcf/current/gcc/LAPACK/lib" "/soft/libraries/alcf/current/gcc/BLAS/lib" "/soft/libraries/petsc/3.5.2.1/xl-opt/lib/" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "metis" "scalapack" "lapack" "blas" "xlf90" "xlopt" "xl" "xlfmath" "xlsmp" "gfortran" "m" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS_CLANG "/soft/libraries/petsc/3.5.2.1/xl-opt/include/" PARENT_SCOPE)
set(QUICC_LIBDIR_MUMPS_CLANG "/soft/compilers/bgclang/xlsmp-nonconflicting/ibmcmp-feb2015" "/soft/compilers/ibmcmp-feb2015/xlf/bg/14.1/bglib64" "/soft/libraries/alcf/current/gcc/LAPACK/lib" "/soft/libraries/alcf/current/gcc/BLAS/lib" "/soft/libraries/alcf/current/gcc/SCALAPACK/lib" "/soft/libraries/petsc/3.5.2.1/xl-opt/lib/" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)
set(QUICC_INCLUDES_HDF5_CLANG "/soft/libraries/hdf5/current/cnk-gcc/current/include/" PARENT_SCOPE)
set(QUICC_LIBDIR_HDF5_CLANG "/soft/libraries/alcf/current/gcc/ZLIB/lib/" "/soft/libraries/hdf5/current/cnk-gcc/current/lib/" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(QUICC_MPIIMPLS "MPICH" PARENT_SCOPE)

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

set(QUICC_CC_SERIAL_CLANG "bgclang++" PARENT_SCOPE)

set(QUICC_CC_MPI_CLANG "mpicxx" PARENT_SCOPE)

set(QUICC_CC_ARCH_CLANG "-O2" PARENT_SCOPE)

set(QUICC_CC_INC_CLANG "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_CLANG ${QUICC_CC_INC_CLANG} PARENT_SCOPE)

set(QUICC_CC_LIB_CLANG "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_CLANG ${QUICC_CC_LIB_CLANG} PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON27 "pthread" "z" "dl" "util" "m" "ssl" "crypto" "python2.7" PARENT_SCOPE)

set(QUICC_LIBDIR_PYTHON27_CLANG "/soft/interpreters/python-2.7.9/powerpc64-bgq-linux/lib/" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27_CLANG "/soft/interpreters/python-2.7.9/powerpc64-bgq-linux/include/python2.7" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES  "" PARENT_SCOPE)
