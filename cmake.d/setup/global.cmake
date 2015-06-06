###################################################
#--------------------- COMPILER ------------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_COMPILERS "Compiler" GEOMHDISCC_COMPILER compilerTest)


###################################################
#----------------- SMART POINTER -----------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_SMARTPTRS "Shared pointer" GEOMHDISCC_SMARTPTR sharedptrTest)

if(sharedptrTest)
   geomhdiscc_add_definition(GEOMHDISCC_SMARTPTR)
endif(sharedptrTest)


###################################################
#------------ THREADS PARALLELISATION ------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_THREADSMODELS "Threads paralellization" GEOMHDISCC_THREADS threadsTest)

if(threadsTest)
   geomhdiscc_add_definition(GEOMHDISCC_THREADS)
endif(threadsTest)


###################################################
#------------- FFT IMPLEMENTATION ----------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_FFTS "FFT implementation" GEOMHDISCC_FFT fftTest)

if(fftTest)
   geomhdiscc_add_definition(GEOMHDISCC_FFT)
endif(fftTest)


###################################################
#------------- FFT PLAN COMPUTATION --------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_FFTPLANS "FFT plan" GEOMHDISCC_FFTPLAN fftplanTest)

if(fftplanTest)
   geomhdiscc_add_definition(GEOMHDISCC_FFTPLAN)
endif(fftplanTest)


###################################################
#--------- LINEAR ALGEBRA IMPLEMENTATION ---------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_LINALGS "Linear algebra" GEOMHDISCC_LINALG linalgTest)

if(linalgTest)
   geomhdiscc_add_definition(GEOMHDISCC_LINALG)
endif(linalgTest)


###################################################
#----- SPARSE LINEAR ALGEBRA IMPLEMENTATION ------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_SPLINALGS "Sparse linear algebra" GEOMHDISCC_SPLINALG splinalgTest)

if(splinalgTest)
   geomhdiscc_add_definition(GEOMHDISCC_SPLINALG)

   if(GEOMHDISCC_SPLINALG STREQUAL "MUMPS")
      option(GEOMHDISCC_MPISPSOLVE "Use MPI sparse solver?" OFF)
      if(GEOMHDISCC_MPISPSOLVE)
         add_definitions("-DGEOMHDISCC_MPISPSOLVE")
      endif(GEOMHDISCC_MPISPSOLVE)
   endif(GEOMHDISCC_SPLINALG STREQUAL "MUMPS")
endif(splinalgTest)


###################################################
#------------------ LARGE IO FORMAT --------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_LARGEIOS "Large IO format" GEOMHDISCC_LARGEIO largeioTest)

if(largeioTest)
   geomhdiscc_add_definition(GEOMHDISCC_LARGEIO)
endif(largeioTest)


###################################################
#--------------- MPI IMPLEMENTATION --------------#
###################################################

geomhdiscc_provide_choice(GEOMHDISCC_MPIIMPLS "MPI implementation" GEOMHDISCC_MPIIMPL mpiimplTest)

if(mpiimplTest)
   geomhdiscc_add_definition(GEOMHDISCC_MPIIMPL)
endif(mpiimplTest)


###################################################
#----------------- SH NORMALIZATION --------------#
###################################################

set(GEOMHDISCC_SHNORM "Unity")
geomhdiscc_add_definition(GEOMHDISCC_SHNORM)


###################################################
#-------------- MULTIPLE PRECISION ---------------#
###################################################

#
# Use multiple precision computation for inititialisation?.
#
option(GEOMHDISCC_MULTPRECISION "Use multiple precision computation for initialisation?" OFF)
if(GEOMHDISCC_MULTPRECISION)
   geomhdiscc_provide_choice(GEOMHDISCC_MPLIBS "Multiple precision" GEOMHDISCC_MPLIB mpTest)
   if(mpTest)
      if(NOT GEOMHDISCC_LINALG STREQUAL "Eigen")
         message(SEND_ERROR "------->>> Can't use multiple precision computations with selected implementation <<<-------")
      else(NOT GEOMHDISCC_LINALG STREQUAL "Eigen")
         add_definitions("-DGEOMHDISCC_MULTPRECISION")
         message(STATUS " --> Multiple precision computations are enabled")
      endif(NOT GEOMHDISCC_LINALG STREQUAL "Eigen")
   endif(mpTest)
endif(GEOMHDISCC_MULTPRECISION)
