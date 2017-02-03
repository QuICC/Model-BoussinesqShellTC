###################################################
#----------- Desktop CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "GCC" "Scalasca" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(QUICC_SMARTPTRS "Boost" "TR1" "cxx0x" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(QUICC_THREADSMODELS "None" "Pthreads" "OpenMP" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(QUICC_FFTS "FFTW" "cuFFT" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW_PTHREADS "fftw3_threads" "fftw3" PARENT_SCOPE)
set(QUICC_LIBRARIES_FFTW_OPENMP "fftw3_omp" "fftw3" PARENT_SCOPE)
set(QUICC_LIBRARIES_CUFFT "cudart" "cufft" "cufftw" PARENT_SCOPE)
set(QUICC_INCLUDES_CUFFT "/opt/cuda/include" "/opt/cuda/sdk/common/inc/" PARENT_SCOPE)
set(QUICC_LIBDIR_CUFFT "/opt/cuda/lib64" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "SparseLU" "MUMPS" "UmfPack" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "MUMPS" "UmfPack" PARENT_SCOPE)
set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "metis" "openblas" "goto2" "mpi" "gfortran" "mpi_mpifh" "cmumps" "smumps" "pord" "scalapack" "parmetis" PARENT_SCOPE)
#set(QUICC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "parmetis" "ptesmumps" "ptscotch" "ptscotcherr" "ptscotchparmetis" "ptscotcherrexit" "scotch" "scotcherr" "scotcherrexit" "pord" "scalapack" "metis" "atllapack" "f77blas" "mpi" "gfortran" "mpi_mpifh" PARENT_SCOPE)
set(QUICC_INCLUDES_MUMPS "/cluster/apps/mumps/5.0.1/x86_64/gcc_4.8.2/openmpi_1.6.5/include" PARENT_SCOPE)
set(QUICC_LIBDIR_MUMPS "/cluster/apps/mumps/5.0.1/x86_64/gcc_4.8.2/openmpi_1.6.5/lib" "/cluster/apps/openmpi/1.6.5/x86_64/gcc_4.8.2/lib" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)

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

set(QUICC_CC_ARCH_GCC "-march=native -O2 -fpermissive" PARENT_SCOPE)

set(QUICC_CC_OPENMP_GCC "-fopenmp" PARENT_SCOPE)

set(QUICC_CC_INC_GCC "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_GCC "${QUICC_CC_INC_GCC}" PARENT_SCOPE)

set(QUICC_CC_LIB_GCC "-rpath" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_GCC "${QUICC_CC_LIB_GCC}" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" "python34" PARENT_SCOPE)

set(QUICC_LIBRARIES_PYTHON27  "$ENV{PYTHON_ROOT}/lib64/libpython2.7.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27  "$ENV{PYTHON_ROOT}/include/python2.7" PARENT_SCOPE)


set(QUICC_LIBRARIES_PYTHON34 "/usr/lib64/libpython3.4.so" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON34 "/usr/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
