###################################################
#--------------------- COMPILER ------------------#
###################################################

quicc_provide_choice(QUICC_COMPILERS "Compiler" QUICC_COMPILER compilerTest)


###################################################
#----------------- SMART POINTER -----------------#
###################################################

if(QUICC_SMARTPTR STREQUAL "")
   list(GET QUICC_SMARTPTRS 0 QUICC_SMARTPTR)
endif(QUICC_SMARTPTR STREQUAL "")
quicc_provide_choice(QUICC_SMARTPTRS "Shared pointer" QUICC_SMARTPTR sharedptrTest)

if(sharedptrTest)
   quicc_add_definition(QUICC_SMARTPTR)
endif(sharedptrTest)


###################################################
#------------ THREADS PARALLELISATION ------------#
###################################################

if(QUICC_THREADS STREQUAL "")
   list(GET QUICC_THREADSMODELS 0 QUICC_THREADS)
endif(QUICC_THREADS STREQUAL "")
quicc_provide_choice(QUICC_THREADSMODELS "Threads paralellization" QUICC_THREADS threadsTest)

if(threadsTest)
   quicc_add_definition(QUICC_THREADS)
endif(threadsTest)


###################################################
#------------- FFT IMPLEMENTATION ----------------#
###################################################

if(QUICC_FFT STREQUAL "")
   list(GET QUICC_FFTS 0 QUICC_FFT)
endif(QUICC_FFT STREQUAL "")
quicc_provide_choice(QUICC_FFTS "FFT implementation" QUICC_FFT fftTest)

if(fftTest)
   quicc_add_definition(QUICC_FFT)
endif(fftTest)


###################################################
#------------- FFT PLAN COMPUTATION --------------#
###################################################
# Default FFT plan is medium
if(QUICC_FFTPLAN STREQUAL "")
  list(GET QUICC_FFTPLANS 0 QUICC_FFTPLAN)
endif(QUICC_FFTPLAN STREQUAL "")
quicc_provide_choice(QUICC_FFTPLANS "FFT plan" QUICC_FFTPLAN fftplanTest)

if(fftplanTest)
   quicc_add_definition(QUICC_FFTPLAN)
endif(fftplanTest)


###################################################
#--------- LINEAR ALGEBRA IMPLEMENTATION ---------#
###################################################

if(QUICC_LINALG STREQUAL "")
   list(GET QUICC_LINALGS 0 QUICC_LINALG)
endif(QUICC_LINALG STREQUAL "")
quicc_provide_choice(QUICC_LINALGS "Linear algebra" QUICC_LINALG linalgTest)

if(linalgTest)
   quicc_add_definition(QUICC_LINALG)
endif(linalgTest)


###################################################
#----- SPARSE LINEAR ALGEBRA IMPLEMENTATION ------#
###################################################

if(QUICC_SPLINALG STREQUAL "")
  list(GET QUICC_SPLINALGS 0 QUICC_SPLINALG)
endif(QUICC_SPLINALG STREQUAL "")
quicc_provide_choice(QUICC_SPLINALGS "Sparse linear algebra" QUICC_SPLINALG splinalgTest)

if(splinalgTest)
   quicc_add_definition(QUICC_SPLINALG)

   if(QUICC_SPLINALG STREQUAL "MUMPS")
      option(QUICC_MPISPSOLVE "Use MPI sparse solver?" OFF)
      if(QUICC_MPISPSOLVE)
         add_definitions("-DQUICC_MPISPSOLVE")
      endif(QUICC_MPISPSOLVE)
   endif(QUICC_SPLINALG STREQUAL "MUMPS")
endif(splinalgTest)


###################################################
#---- SPARSE SPD LINEAR ALGEBRA IMPLEMENTATION ---#
###################################################

if(QUICC_SPSPDLINALG STREQUAL "")
   set(QUICC_SPSPDLINALG "SimplicialLDLT")
endif(QUICC_SPSPDLINALG STREQUAL "")
mark_as_advanced(FORCE QUICC_SPSPDLINALG)
quicc_provide_choice(QUICC_SPSPDLINALGS "Sparse SPD linear algebra" QUICC_SPSPDLINALG spspdlinalgTest)

if(spspdlinalgTest)
   quicc_add_definition(QUICC_SPSPDLINALG)
endif(spspdlinalgTest)


###################################################
#---- SPARSE TRI LINEAR ALGEBRA IMPLEMENTATION ---#
###################################################

if(QUICC_SPTRILINALG STREQUAL "")
   set(QUICC_SPTRILINALG "SparseLU")
endif(QUICC_SPTRILINALG STREQUAL "")
mark_as_advanced(FORCE QUICC_SPTRILINALG)
quicc_provide_choice(QUICC_SPTRILINALGS "Sparse triangular linear algebra" QUICC_SPTRILINALG sptrilinalgTest)

if(sptrilinalgTest)
   quicc_add_definition(QUICC_SPTRILINALG)
endif(sptrilinalgTest)


###################################################
#------------------ LARGE IO FORMAT --------------#
###################################################

if(QUICC_LARGEIOS STREQUAL "")
   list(GET QUICC_LARGEIOS 0 QUICC_LARGEIO)
endif(QUICC_LARGEIOS STREQUAL "")
quicc_provide_choice(QUICC_LARGEIOS "Large IO format" QUICC_LARGEIO largeioTest)

if(largeioTest)
   quicc_add_definition(QUICC_LARGEIO)
endif(largeioTest)


###################################################
#-------------- HDF5 COMPLEX FORMAT --------------#
###################################################

set(QUICC_HDF5_CMPLXS "Array" "Struct")
if(QUICC_HDF5_CMPLX STREQUAL "")
   set(QUICC_HDF5_CMPLX "Array")
endif(QUICC_HDF5_CMPLX STREQUAL "")
mark_as_advanced(FORCE QUICC_HDF5_CMPLX)
quicc_provide_choice(QUICC_HDF5_CMPLXS "HDF5 complex format" QUICC_HDF5_CMPLX h5cmplxTest)

if(h5cmplxTest)
   quicc_add_definition(QUICC_HDF5_CMPLX)
endif(h5cmplxTest)


###################################################
#--------------- MPI IMPLEMENTATION --------------#
###################################################

quicc_provide_choice(QUICC_MPIIMPLS "MPI implementation" QUICC_MPIIMPL mpiimplTest)

if(mpiimplTest)
   quicc_add_definition(QUICC_MPIIMPL)
endif(mpiimplTest)


###################################################
#-------- ASSOCIATED LEGENDRE TRANSFORM ----------#
###################################################

set(QUICC_ALEGTRAS "Matrix" "Fly")
if(QUICC_ALEGTRA STREQUAL "")
   set(QUICC_ALEGTRA "Matrix")
endif(QUICC_ALEGTRA STREQUAL "")
mark_as_advanced(FORCE QUICC_ALEGTRA)
quicc_provide_choice(QUICC_ALEGTRAS "Associated Legendre implementation" QUICC_ALEGTRA alegtraTest)

if(alegtraTest)
   quicc_add_definition(QUICC_ALEGTRA)
endif(alegtraTest)


###################################################
#----------------- SH NORMALIZATION --------------#
###################################################

set(QUICC_SHNORMS "Unity" "Schmidt")
if(QUICC_SHNORM STREQUAL "")
   set(QUICC_SHNORM "Unity")
endif(QUICC_SHNORM STREQUAL "")
mark_as_advanced(FORCE QUICC_SHNORM)
quicc_provide_choice(QUICC_SHNORMS "Spherical harmonics normalization" QUICC_SHNORM shnormTest)

if(shnormTest)
   quicc_add_definition(QUICC_SHNORM)
endif(shnormTest)


###################################################
#-------------- WORLAND TRANSFORM ----------------#
###################################################

set(QUICC_WORLANDTRAS "Matrix" "Fly")
if(QUICC_WORLANDTRA STREQUAL "")
   set(QUICC_WORLANDTRA "Matrix")
endif(QUICC_WORLANDTRA STREQUAL "")
mark_as_advanced(FORCE QUICC_WORLANDTRA)
quicc_provide_choice(QUICC_WORLANDTRAS "Worland implementation" QUICC_WORLANDTRA worlandtraTest)

if(worlandtraTest)
   quicc_add_definition(QUICC_WORLANDTRA)
endif(worlandtraTest)


###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

#
# Use multiple precision computation for inititialisation?.
#
option(QUICC_MULTPRECISION "Enable multiple precision computations?" OFF)
if(QUICC_MULTPRECISION)
   quicc_provide_choice(QUICC_MPBACKENDS "Multiple precision backend" QUICC_MPBACKEND mpTest)
   if(mpTest)
      if(NOT QUICC_LINALG STREQUAL "Eigen")
         message(SEND_ERROR "------->>> Can't use multiple precision computations with selected implementation <<<-------")
      else(NOT QUICC_LINALG STREQUAL "Eigen")
         add_definitions("-DQUICC_MULTPRECISION")
         quicc_add_definition(QUICC_MPBACKEND)
         if(NOT QUICC_MPBACKEND STREQUAL "quad")
           if(NOT DEFINED QUICC_MULTPRECISION_DIGITS)
              set(QUICC_MULTPRECISION_DIGITS 50)
           endif(NOT DEFINED QUICC_MULTPRECISION_DIGITS)
           set(QUICC_MULTPRECISION_DIGITS ${QUICC_MULTPRECISION_DIGITS} CACHE STRING "Multiple precision digits" FORCE)
         add_definitions("-DQUICC_MULTPRECISION_DIGITS=${QUICC_MULTPRECISION_DIGITS}")
         endif(NOT QUICC_MPBACKEND STREQUAL "quad")
         message(STATUS " --> Multiple precision computations are enabled")
      endif(NOT QUICC_LINALG STREQUAL "Eigen")
   endif(mpTest)
endif(QUICC_MULTPRECISION)

###################################################
#---------- TRANSFORM TREE OPTIMIZATION ----------#
###################################################

# Disable by default as it doesn't work for all cases
option(QUICC_OPTIMIZE_TREE "Optimize transform tree?" ON)
mark_as_advanced(FORCE QUICC_OPTIMIZE_TREE)
if(QUICC_OPTIMIZE_TREE)
   add_definitions("-DQUICC_OPTIMIZE_TREE")
endif(QUICC_OPTIMIZE_TREE)


###################################################
#------------------ EMBEDDED PYTHON --------------#
###################################################

quicc_provide_choice(QUICC_PYTHONS "Python version" QUICC_PYTHON pythonTest)
