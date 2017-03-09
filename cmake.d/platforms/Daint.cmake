
###################################################
#----------- Pleiades CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "GCC" PARENT_SCOPE)
#set(QUICC_DISABLE_RDYNAMIC ON PARENT_SCOPE)
set(QUICC_ENABLE_DYNAMIC ON PARENT_SCOPE)

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
set(QUICC_INCLUDES_FFTW "$ENV{FFTW_INC}" PARENT_SCOPE)
set(QUICC_LIBDIR_FFTW "$ENV{FFTW_DIR}" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "SparseLU" "MUMPS" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS "$ENV{CRAY_TPSL_PREFIX_DIR}/include" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "MUMPS" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "hdf5" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(QUICC_MPIIMPLS "CRAYMPICH" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(QUICC_MPLIBS "mpfr" PARENT_SCOPE)
set(QUICC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(QUICC_CC_SERIAL_GCC "CC" PARENT_SCOPE)

set(QUICC_CC_MPI_GCC "CC" PARENT_SCOPE)

set(QUICC_CC_ARCH_GCC "-O2" PARENT_SCOPE)

set(QUICC_CC_INC_GCC "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_GCC ${QUICC_CC_INC_GCC} PARENT_SCOPE)

set(QUICC_CC_LIB_GCC "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_GCC "${QUICC_CC_LIB_GCC}" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python35" PARENT_SCOPE)


set(QUICC_LIBRARIES_PYTHON27 "python2.7" "util" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27 "$ENV{EBROOTPYTHON}/include/python2.7" PARENT_SCOPE)
set(QUICC_LIBDIR_PYTHON27 "$ENV{EBROOTPYTHON}/lib" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON34 "python3.4m" "util" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON34 "$ENV{EBROOTPYTHON}/include/python3.4m" PARENT_SCOPE)
set(QUICC_LIBDIR_PYTHON34 "$ENV{EBROOTPYTHON}/lib" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON35 "python3.5m" "util" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON35 "$ENV{EBROOTPYTHON}/include/python3.5m" "$ENV{EBROOTPYTHON}/lib/python3.5/site-packages/numpy-1.11.2-py3.5-linux-x86_64.egg/numpy/core/include"  PARENT_SCOPE)
set(QUICC_LIBDIR_PYTHON35 "$ENV{EBROOTPYTHON}/lib" PARENT_SCOPE)


###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES  "" PARENT_SCOPE)
