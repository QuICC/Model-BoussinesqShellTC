###################################################
#----------- Janus CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(GEOMHDISCC_COMPILERS "GCC" "Intel" PARENT_SCOPE)

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
set(GEOMHDISCC_INCLUDES_CUFFT "/curc/tools/x_86_64/rh6/cuda/5.5/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_CUFFT "/curc/tools/x_86_64/rh6/cuda/5.5/lib64" PARENT_SCOPE)

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

set(GEOMHDISCC_SPLINALGS "SuperLU" "UmfPack" "MUMPS" "SparseLU" "KentLU" "MKLPardiso" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SUPERLU "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "superlu_4.3" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MUMPS "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "dmumps" "zmumps" "mumps_common" "pord" "scalapack" "ptesmumps" "esmumps" "ptscotch" "scotch" "scotcherr" "scotcherrexit" "parmetis" "metis" "ifcore" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_KENTLU "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "klu" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MKLPARDISO "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_UMFPACK_GCC "/home/phma6156/share/gcc/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_UMFPACK_GCC "/home/phma6156/share/gcc/lib" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SUPERLU_GCC "/home/phma6156/share/gcc/SuperLU_4.3/SRC" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_SUPERLU_GCC "/home/phma6156/share/gcc/SuperLU_4.3/lib" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_UMFPACK_INTEL "/home/phma6156/share/intel/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_UMFPACK_INTEL "/home/phma6156/share/intel/lib" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/SRC" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/lib" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_MUMPS_INTEL "/home/phma6156/share/intel/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_MUMPS_INTEL "/home/phma6156/share/intel/lib" PARENT_SCOPE)

###################################################
#--------- AVAILABLE TRANSFORM OPERATORS ---------#
###################################################

set(GEOMHDISCC_TRANSOPS "Forward" "Recurrence" "Backward" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)

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

set(GEOMHDISCC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_INTEL "mpic++" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_INTEL "-O2" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_INTEL ${GEOMHDISCC_CC_INC_INTEL} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_INTEL ${GEOMHDISCC_CC_LIB_INTEL} PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

# Find python headers and library
set(PYTHON_LIBRARY "/curc/tools/x_86_64/rh6/anaconda/2.0.0/lib/libpython2.7.so")
set(PYTHON_INCLUDE_DIR "/curc/tools/x_86_64/rh6/anaconda/2.0.0/include/python2.7")

set(GEOMHDISCC_LIBRARIES ${PYTHON_LIBRARY} PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES ${PYTHON_INCLUDE_DIR} PARENT_SCOPE)
