###################################################
#----------- Pleiades CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(GEOMHDISCC_COMPILERS "Intel" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(GEOMHDISCC_SMARTPTRS "TR1" "cxx0x" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(GEOMHDISCC_THREADSMODELS "None" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(GEOMHDISCC_FFTS "FFTW" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_FFTW "/nobackupp9/nfeather/marti/libs/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_FFTW "/nobackupp9/nfeather/marti/libs/lib" PARENT_SCOPE)

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

set(GEOMHDISCC_SPLINALGS "MUMPS" "SparseLU" "UmfPack" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_MUMPS "/nobackupp9/nfeather/marti/libs/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_MUMPS "/nobackupp9/nfeather/marti/libs/lib" "$ENV{MKLROOT}/lib/intel64" PARENT_SCOPE)
#set(GEOMHDISCC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "mkl_scalapack_lp64" "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "mpifort" "ifcore" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "mkl_scalapack_lp64" "mkl_intel_lp64" "mkl_core" "mkl_sequential" "mkl_blacs_sgimpt_lp64" "pthread" "m" "dl" "ifcore" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
#set(GEOMHDISCC_LIBRARIES_SUPERLU "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "superlu_4.3" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_UMFPACK_INTEL "$ENV{HOME}/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_UMFPACK_INTEL $ENV{TACC_MKL_LIB} "$ENV{HOME}/local/lib" PARENT_SCOPE)
#set(GEOMHDISCC_INCLUDES_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/SRC" PARENT_SCOPE)
#set(GEOMHDISCC_LIBDIR_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/lib" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)
#set(GEOMHDISCC_INCLUDES_HDF5 "" PARENT_SCOPE)
#set(GEOMHDISCC_LIBDIR_HDF5 "" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(GEOMHDISCC_MPIIMPLS "SGIMPT" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(GEOMHDISCC_MPLIBS "mpfr" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(GEOMHDISCC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_INTEL "icpc" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_INTEL "-O2 -xAVX" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_INTEL ${GEOMHDISCC_CC_INC_INTEL} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_INTEL "${GEOMHDISCC_CC_LIB_INTEL}" "mpi" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_PYTHONS "python34" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON34 "pthread" "dl" "util" "m" "python3.4" PARENT_SCOPE)
#set(GEOMHDISCC_LIBDIR_PYTHON27 "$ENV{INCLUDE}" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON34 "/nasa/python/3.4.3/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_LIBRARIES "" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES  "" PARENT_SCOPE)
