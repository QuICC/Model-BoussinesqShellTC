###################################################
#----------- Desktop CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(GEOMHDISCC_COMPILERS "GCC" "Scalasca" PARENT_SCOPE)

###################################################
#----------- AVAILABLE SMART POINTERS ------------#
###################################################

set(GEOMHDISCC_SMARTPTRS "Boost" "TR1" "cxx0x" PARENT_SCOPE)

###################################################
#----------- AVAILABLE THREADS MODELS ------------#
###################################################

set(GEOMHDISCC_THREADSMODELS "None" "Pthreads" "OpenMP" PARENT_SCOPE)

###################################################
#----------- AVAILABLE FFT VERSIONS --------------#
###################################################

set(GEOMHDISCC_FFTS "FFTW" "cuFFT" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW_PTHREADS "fftw3_threads" "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_FFTW_OPENMP "fftw3_omp" "fftw3" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_CUFFT "cudart" "cufft" "cufftw" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_CUFFT "/opt/cuda/include" "/opt/cuda/sdk/common/inc/" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_CUFFT "/opt/cuda/lib64" PARENT_SCOPE)

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

set(GEOMHDISCC_SPLINALGS "SuperLU" "UmfPack" "SparseLU" "MUMPS" "Pardiso" "SPQR" "SparseQR" "BiCGSTAB" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SUPERLU "superlu" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SUPERLU "/usr/include/superlu" PARENT_SCOPE)
#set(GEOMHDISCC_LIBRARIES_SPARSELU "/usr/local/lib/libmetis.a" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_SPARSELU "/usr/local/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MUMPS "dmumps" "zmumps" "mumps_common" "parmetis" "ptesmumps" "ptscotch" "ptscotcherr" "ptscotchparmetis" "ptscotcherrexit" "scotch" "scotcherr" "scotcherrexit" "pord" "scalapack" "metis" "atllapack" "f77blas" "mpi" "gfortran" "mpi_mpifh" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_MUMPS "/usr/local/share/petsc-3.6.3_complex/arch-linux2-c-opt/include" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_MUMPS "/usr/local/share/petsc-3.6.3_complex/arch-linux2-c-opt/lib" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_PARDISO "pardiso412-GNU450-X86-64" PARENT_SCOPE)
set(GEOMHDISCC_LIBDIR_PARDISO "/usr/local/share/Pardiso/lib" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(GEOMHDISCC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "UmfPack" "SparseLU" "MUMPS" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(GEOMHDISCC_SPTRILINALGS "SparseLU" "UmfPack" "MUMPS" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(GEOMHDISCC_LARGEIOS "HDF5" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)

###################################################
#-------------- MPI IMPLEMENTATION ---------------#
###################################################

set(GEOMHDISCC_MPIIMPLS "OpenMPI" PARENT_SCOPE)

###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

set(GEOMHDISCC_MPLIBS "mpfr" PARENT_SCOPE)
set(GEOMHDISCC_LIBRARIES_MPFR "mpfr" PARENT_SCOPE)

###################################################
#--------------- COMPILER SETTINGS ---------------#
###################################################

set(GEOMHDISCC_CC_SERIAL_GCC "g++" PARENT_SCOPE)
set(GEOMHDISCC_CC_SERIAL_SCALASCA "g++" PARENT_SCOPE)

set(GEOMHDISCC_CC_MPI_GCC "mpic++" PARENT_SCOPE)
set(GEOMHDISCC_CC_MPI_SCALASCA "mpic++" PARENT_SCOPE)

set(GEOMHDISCC_CC_ARCH_GCC "-march=native -O2" PARENT_SCOPE)
set(GEOMHDISCC_CC_ARCH_SCALASCA "-march=native -O2" PARENT_SCOPE)

set(GEOMHDISCC_CC_OPENMP_GCC "-fopenmp" PARENT_SCOPE)
set(GEOMHDISCC_CC_OPENMP_SCALASCA "-fopenmp" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_SCALASCA "" PARENT_SCOPE)

set(GEOMHDISCC_CC_INC_MPI_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_INC_MPI_SCALASCA "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_SCALASCA "" PARENT_SCOPE)

set(GEOMHDISCC_CC_LIB_MPI_GCC "" PARENT_SCOPE)
set(GEOMHDISCC_CC_LIB_MPI_SCALASCA "" PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_PYTHONS "python27" "python34" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON27 "/usr/lib64/libpython2.7.so" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON27 "/usr/include/python2.7" PARENT_SCOPE)

set(GEOMHDISCC_LIBRARIES_PYTHON34 "/usr/lib64/libpython3.4.so" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES_PYTHON34 "/usr/include/python3.4" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(GEOMHDISCC_LIBRARIES "" PARENT_SCOPE)
set(GEOMHDISCC_INCLUDES "" PARENT_SCOPE)
