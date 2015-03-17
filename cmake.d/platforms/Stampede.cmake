###################################################
#----------- Janus CMAKE SETUP FILE ------------#
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
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(GEOMHDISCC_FFTS "FFTW" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_FFTW $ENV{TACC_FFTW3_INC} PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_FFTW $ENV{TACC_FFTW3_LIB} PARENT_SCOPE)

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
#set(GEOMHDISCC_INCLUDES_MUMPS $ENV{TACC_MUMPS_INC} PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_MUMPS $ENV{TACC_MUMPS_INC} PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_MUMPS $ENV{TACC_MUMPS_LIB} $ENV{TACC_MKL_LIB} PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "pord" "parmetis" "scalapack" "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "umfpack" "amd" "cholmod" "ccolamd" "colamd" "camd" "metis" "suitesparseconfig" PARENT_SCOPE)
#set(GEOMHDISCC_LIBRARIES_SUPERLU "mkl_intel_lp64" "mkl_sequential" "mkl_core" "pthread" "m" "superlu_4.3" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_UMFPACK_INTEL "$ENV{HOME}/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_UMFPACK_INTEL $ENV{TACC_MKL_LIB} "$ENV{HOME}/local/lib" PARENT_SCOPE)
#set(GEOMHDISCC_INCLUDES_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/SRC" PARENT_SCOPE)
#set(GEOMHDISCC_LIBDIR_SUPERLU_INTEL "/home/phma6156/share/intel/SuperLU_4.3/lib" PARENT_SCOPE)

###################################################
#--------- AVAILABLE TRANSFORM OPERATORS ---------#
###################################################

set(GEOMHDISCC_TRANSOPS "Forward" "Recurrence" "Backward" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "z" "hdf5" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_HDF5 $ENV{TACC_HDF5_INC} PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_HDF5 "$ENV{TACC_HDF5_LIB}" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(GEOMHDISCC_MPIIMPLS "MVAPICH" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(GEOMHDISCC_MPLIBS "mpfr" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(GEOMHDISCC_CC_SERIAL_INTEL "icpc" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_INTEL "mpicxx" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_INTEL "-xhost -O2" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_INTEL ${GEOMHDISCC_CC_INC_INTEL} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_INTEL ${GEOMHDISCC_CC_LIB_INTEL} PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

# Find python headers and library
set(PYTHON_LIBRARY "pthread" "dl" "util" "m" "$ENV{TACC_PYTHON_LIB}/libpython2.7.so")
set(PYTHON_INCLUDE_DIR "$ENV{TACC_PYTHON_DIR}/include/python2.7")

set(GEOMHDISCC_LIBRARIES ${PYTHON_LIBRARY} PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES ${PYTHON_INCLUDE_DIR} PARENT_SCOPE)
