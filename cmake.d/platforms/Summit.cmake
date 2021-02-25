###################################################
#----------- Janus CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "Intel" PARENT_SCOPE)

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
set(QUICC_INCLUDES_FFTW "$ENV{CURC_FFTW_INC}" PARENT_SCOPE)
set(QUICC_LIBDIR_FFTW "$ENV{CURC_FFTW_LIB}" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "UmfPack" "MUMPS"  "SparseLU" PARENT_SCOPE)
set(QUICC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "mkl_scalapack_ilp64" "mkl_intel_lp64" "mkl_sequential" "mkl_core" "mkl_blacs_intelmpi_ilp64" "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "metis" "ifcore" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_GCC "/home/mica5951/SuiteSparse/include" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_GCC "/home/mica5951/SuiteSparse/lib" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_INTEL "/home/mica5951/SuiteSparse/include" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_INTEL "/home/mica5951/SuiteSparse/lib" "$ENV{CURC_MKL_LIB}" PARENT_SCOPE)
#set(QUICC_INCLUDES_MUMPS_INTEL "/home/phma6156/share/intel/mumps_4/include" PARENT_SCOPE)
#set(QUICC_LIBDIR_MUMPS_INTEL "/home/phma6156/share/intel/mumps_4/lib" "/home/phma6156/share/intel/scalapack/lib" "/home/phma6156/share/intel/scotch_5/lib" "/home/phma6156/share/intel/parmetis_3/lib" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS "$ENV{CURC_PETSC_INC}" PARENT_SCOPE)
set(QUICC_LIBDIR_MUMPS "$ENV{CURC_PETSC_LIB}" "$ENV{CURC_MKL_LIB}" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "UmfPack" "MUMPS" "SparseLU" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)
set(QUICC_INCLUDES_HDF5 "$ENV{CURC_HDF5_INC}" PARENT_SCOPE)
set(QUICC_LIBDIR_HDF5 "$ENV{CURC_HDF5_LIB}" PARENT_SCOPE)

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

set(QUICC_CC_ARCH_GCC "-std=c++11 -march=native -O2" PARENT_SCOPE)

set(QUICC_CC_INC_GCC "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_GCC ${QUICC_CC_INC_GCC} PARENT_SCOPE)

set(QUICC_CC_LIB_GCC "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_GCC ${QUICC_CC_LIB_GCC} PARENT_SCOPE)

set(QUICC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(QUICC_CC_MPI_INTEL "mpicxx" PARENT_SCOPE)

set(QUICC_CC_ARCH_INTEL "-std=c++11 -O2 -xHost" PARENT_SCOPE)

set(QUICC_CC_INC_INTEL "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_INTEL ${QUICC_CC_INC_INTEL} PARENT_SCOPE)

set(QUICC_CC_LIB_INTEL "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_INTEL ${QUICC_CC_LIB_INTEL} PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" "python34" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON27 "$ENV{CURC_PYTHON_LIB}/libpython2.7.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27 "$ENV{CURC_PYTHON_INC}/python2.7" PARENT_SCOPE)
set(QUICC_LIBRARIES_PYTHON34 "util" "$ENV{CURC_PYTHON_LIB}/libpython3.4m.a" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON34 "$ENV{CURC_PYTHON_INC}/python3.4m" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
