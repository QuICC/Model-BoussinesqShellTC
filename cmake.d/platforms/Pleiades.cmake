###################################################
#----------- Pleiades CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "Intel" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(QUICC_SMARTPTRS "TR1" "cxx0x" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(QUICC_THREADSMODELS "None" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(QUICC_FFTS "FFTW" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(QUICC_INCLUDES_FFTW "/nobackup/mplumle1/libs/include" PARENT_SCOPE)
set(QUICC_LIBDIR_FFTW "/nobackup/mplumle1/libs/lib" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "MUMPS" "SparseLU" "UmfPack" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS "/nobackup/mplumle1/libs/include" PARENT_SCOPE)
set(QUICC_LIBDIR_MUMPS "/nobackup/mplumle1/libs/lib" "$ENV{MKLROOT}/lib/intel64" PARENT_SCOPE)
#set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "mkl_scalapack_lp64" "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "mpifort" "ifcore" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "metis" "mkl_scalapack_lp64" "mkl_intel_lp64" "mkl_core" "mkl_sequential" "mkl_blacs_sgimpt_lp64" "pthread" "m" "dl" "ifcore" PARENT_SCOPE)
set(QUICC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_INTEL "/nobackup/mplumle1/libs/include" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_INTEL "$ENV{MKLROOT}/lib/intel64" "/nobackup/mplumle1/libs/lib" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "MUMPS" "SparseLU" "UmfPack" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "MUMPS" "UmfPack" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)
#set(QUICC_INCLUDES_HDF5 "" PARENT_SCOPE)
#set(QUICC_LIBDIR_HDF5 "" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(QUICC_MPIIMPLS "SGIMPT" PARENT_SCOPE)

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

set(QUICC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(QUICC_CC_MPI_INTEL "icpc" PARENT_SCOPE)

set(QUICC_CC_ARCH_INTEL "-O2 -xAVX" PARENT_SCOPE)

set(QUICC_CC_INC_INTEL "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_INTEL ${QUICC_CC_INC_INTEL} PARENT_SCOPE)

set(QUICC_CC_LIB_INTEL "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_INTEL "${QUICC_CC_LIB_INTEL}" "mpi" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python34" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON34 "pthread" "dl" "util" "m" "python3.4" PARENT_SCOPE)
#set(QUICC_LIBDIR_PYTHON27 "$ENV{INCLUDE}" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON34 "/nasa/python/3.4.3/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES  "" PARENT_SCOPE)
