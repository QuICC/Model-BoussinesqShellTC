###################################################
#----------- Janus CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "GCC" "Intel" PARENT_SCOPE)

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

set(QUICC_FFTS "FFTW" "cuFFT" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(QUICC_LIBRARIES_CUFFT "cudart" "cufft" "cufftw" PARENT_SCOPE)
set(QUICC_INCLUDES_CUFFT "/curc/tools/x_86_64/rh6/cuda/5.5/include" PARENT_SCOPE)
set(QUICC_LIBDIR_CUFFT "/curc/tools/x_86_64/rh6/cuda/5.5/lib64" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "UmfPack" "MUMPS"  "SparseLU" "KentLU" "MKLPardiso" PARENT_SCOPE)
set(QUICC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "dmumps" "zmumps" "mumps_common" "pord" "scalapack" "ptesmumps" "esmumps" "ptscotch" "scotch" "scotcherr" "scotcherrexit" "parmetis" "metis" "ifcore" PARENT_SCOPE)
set(QUICC_LIBRARIES_KENTLU "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "klu" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(QUICC_LIBRARIES_MKLPARDISO "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_GCC "/home/phma6156/share/gcc/include" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_GCC "/home/phma6156/share/gcc/lib" PARENT_SCOPE)
set(QUICC_INCLUDES_SUPERLU_GCC "/home/phma6156/share/gcc/SuperLU_4.3/SRC" PARENT_SCOPE)
set(QUICC_LIBDIR_SUPERLU_GCC "/home/phma6156/share/gcc/SuperLU_4.3/lib" PARENT_SCOPE)
set(QUICC_INCLUDES_UMFPACK_INTEL "/home/phma6156/share/intel/suitesparse/include" PARENT_SCOPE)
set(QUICC_LIBDIR_UMFPACK_INTEL "/home/phma6156/share/intel/suitesparse/lib" PARENT_SCOPE)
#set(QUICC_INCLUDES_MUMPS_INTEL "/home/phma6156/share/intel/mumps_4/include" PARENT_SCOPE)
#set(QUICC_LIBDIR_MUMPS_INTEL "/home/phma6156/share/intel/mumps_4/lib" "/home/phma6156/share/intel/scalapack/lib" "/home/phma6156/share/intel/scotch_5/lib" "/home/phma6156/share/intel/parmetis_3/lib" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS_INTEL "/home/phma6156/share/intel/mumps_5/include" PARENT_SCOPE)
set(QUICC_LIBDIR_MUMPS_INTEL "/home/phma6156/share/intel/mumps_5/lib" "/home/phma6156/share/intel/scalapack/lib" "/home/phma6156/share/intel/scotch_6/lib" "/home/phma6156/share/intel/parmetis_4/lib" PARENT_SCOPE)

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

set(QUICC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(QUICC_CC_MPI_INTEL "mpic++" PARENT_SCOPE)

set(QUICC_CC_ARCH_INTEL "-O2 -xHost" PARENT_SCOPE)

set(QUICC_CC_INC_INTEL "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_INTEL ${QUICC_CC_INC_INTEL} PARENT_SCOPE)

set(QUICC_CC_LIB_INTEL "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_INTEL ${QUICC_CC_LIB_INTEL} PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON27 "/curc/tools/x_86_64/rh6/anaconda/2.0.0/lib/libpython2.7.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27 "/curc/tools/x_86_64/rh6/anaconda/2.0.0/include/python2.7" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
