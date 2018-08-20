###################################################
#----------- Desktop CMAKE SETUP FILE ------------#
###################################################

###################################################
#-------------- AVAILABLE COMPILERS --------------#
###################################################

set(QUICC_COMPILERS "CLANG" PARENT_SCOPE)

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

set(QUICC_SPLINALGS "SuperLU" "UmfPack" "SparseLU" PARENT_SCOPE)
set(QUICC_LIBRARIES_UMFPACK "umfpack" PARENT_SCOPE)
set(QUICC_LIBRARIES_SUPERLU "superlu" PARENT_SCOPE)
set(QUICC_INCLUDES_SUPERLU "/usr/include/superlu" PARENT_SCOPE)
set(QUICC_LIBRARIES_SPQR "spqr" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE SPD LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPSPDLINALGS "SimplicialLDLT" "SimplicialLLT" "SparseLU" "SuperLU" "UmfPack" PARENT_SCOPE)

###################################################
#- AVAILABLE SPARSE TRI LINEAR ALGEBRA LIBRARIES -#
###################################################

set(QUICC_SPTRILINALGS "SparseLU" "SuperLU" "UmfPack" PARENT_SCOPE)

###################################################
#----------- AVAILABLE LARGE IO FORMAT -----------#
###################################################

set(QUICC_LARGEIOS "HDF5" PARENT_SCOPE)
#set(QUICC_LIBRARIES_HDF5 "rt" "hdf5" PARENT_SCOPE)
set(QUICC_LIBRARIES_HDF5 "hdf5" PARENT_SCOPE)

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

set(QUICC_CC_SERIAL_CLANG "g++-5" PARENT_SCOPE)

set(QUICC_CC_MPI_CLANG "mpic++" PARENT_SCOPE)

# for debug
set(QUICC_CC_ARCH_GCC "-march=native -O0 -g -stdlib=libstdc++ --std=c++11" PARENT_SCOPE)
# for fast, optimised runs
#set(QUICC_CC_ARCH_GCC "-march=native -O2 -stdlib=libstdc++ --std=c++11" PARENT_SCOPE)

set(QUICC_CC_INC_CLANG "" PARENT_SCOPE)

set(QUICC_CC_INC_MPI_CLANG ${QUICC_CC_INC_CLANG} PARENT_SCOPE)

set(QUICC_CC_LIB_CLANG "" PARENT_SCOPE)

set(QUICC_CC_LIB_MPI_CLANG ${QUICC_CC_LIB_CLANG} PARENT_SCOPE)

###################################################
#--------------- PYTHON LIBRARIES ----------------#
###################################################

set(QUICC_PYTHONS "python27" PARENT_SCOPE)

#set(QUICC_LIBRARIES_PYTHON27 "/usr/lib/libpython2.7.dylib" PARENT_SCOPE)
#set(QUICC_INCLUDES_PYTHON27 "/usr/include/python2.7" PARENT_SCOPE)
set(QUICC_LIBRARIES_PYTHON27 "/Users/stefanomaffei/anaconda/lib/libpython2.7.dylib" PARENT_SCOPE)
set(QUICC_INCLUDES_PYTHON27 "/Users/stefanomaffei/anaconda/include/python2.7" PARENT_SCOPE)

###################################################
#-------------- GENERAL LIBRARIES ----------------#
###################################################

set(QUICC_LIBRARIES "-lboost_math_tr1" PARENT_SCOPE)
set(QUICC_INCLUDES "" PARENT_SCOPE)
