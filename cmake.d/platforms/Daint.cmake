###################################################
#----------- Pleiades CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(GEOMHDISCC_COMPILERS "Cray" "GCC" "Intel" PARENT_SCOPE)
set(GEOMHDISCC_DISABLE_DYNAMIC ON PARENT_SCOPE)

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
set(GEOMHDISCC_INCLUDES_FFTW "$ENV{FFTW_INC}" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_FFTW "$ENV{FFTW_DIR}" PARENT_SCOPE)

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

set(GEOMHDISCC_SPLINALGS "SparseLU" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(GEOMHDISCC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(GEOMHDISCC_SPTRILINALGS "SparseLU" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "hdf5" PARENT_SCOPE)
#set(GEOMHDISCC_INCLUDES_HDF5 "" PARENT_SCOPE)
#set(GEOMHDISCC_LIBDIR_HDF5 "" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(GEOMHDISCC_MPIIMPLS "MPICH2" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(GEOMHDISCC_MPLIBS "mpfr" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(GEOMHDISCC_CC_SERIAL_CRAY "CC" PARENT_SCOPE)
set(GEOMHDISCC_CC_SERIAL_GCC "CC" PARENT_SCOPE)
set(GEOMHDISCC_CC_SERIAL_INTEL "CC" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_CRAY "CC" PARENT_SCOPE)
set(GEOMHDISCC_CC_MPI_GCC "CC" PARENT_SCOPE)
set(GEOMHDISCC_CC_MPI_INTEL "CC" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_CRAY "" PARENT_SCOPE)
set(GEOMHDISCC_CC_ARCH_GCC "-O2" PARENT_SCOPE)
set(GEOMHDISCC_CC_ARCH_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_CRAY "" PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_CRAY ${GEOMHDISCC_CC_INC_CRAY} PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_MPI_GCC ${GEOMHDISCC_CC_INC_GCC} PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_MPI_INTEL ${GEOMHDISCC_CC_INC_INTEL} PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_CRAY "" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_INTEL "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_CRAY "${GEOMHDISCC_CC_LIB_CRAY}" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_MPI_GCC "${GEOMHDISCC_CC_LIB_GCC}" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_MPI_INTEL "${GEOMHDISCC_CC_LIB_INTEL}" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_PYTHONS "python27" "python34" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON27 "python2.7" "util" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON27 "$ENV{PYTHON_DIR}/include/python2.7" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_PYTHON27 "$ENV{PYTHONPATH}" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON34 "python3.4m" "util" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON34 "$ENV{PYTHON_DIR}/include/python3.4m" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_PYTHON34 "$ENV{PYTHONPATH}" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_LIBRARIES "" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES  "" PARENT_SCOPE)
