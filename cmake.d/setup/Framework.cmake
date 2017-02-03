#
# Include function library
#
include(cmake.d/functions.cmake)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

message(STATUS "***********************************************")
message(STATUS "*************** Framework setup ***************")
message(STATUS "***********************************************")

###################################################
#----------------- MEMORY USAGE ------------------#
###################################################

#
# Choose the type of memory usage setup the code is using.
# Possible options are: High, Limited
#
set(QUICC_MEMORYUSAGES "High" "Limited")

quicc_provide_choice(QUICC_MEMORYUSAGES "Memory usage" QUICC_MEMORYUSAGE memoryTest)

if(memoryTest)
   quicc_add_definition(QUICC_MEMORYUSAGE)
endif(memoryTest)


###################################################
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Choose the type of MPI parallelisation or serial setup.
# Possible options are: Serial, Auto, Single1D, Single2D, Tubular
#
set(QUICC_MPIALGOS "Serial" "Auto" "Single1D" "Single2D" "Tubular" "Coupled2D")

quicc_provide_choice(QUICC_MPIALGOS "MPI algorithm" QUICC_MPIALGO mpiTest)

if(mpiTest)
   if(NOT QUICC_MPIALGO STREQUAL "Serial")
      set(QUICC_MPI ON)
      add_definitions("-DQUICC_MPI")
   endif(NOT QUICC_MPIALGO STREQUAL "Serial")

   if(NOT QUICC_MPIALGO STREQUAL "Auto")
      quicc_add_definition(QUICC_MPIALGO)
   else(NOT QUICC_MPIALGO STREQUAL "Auto")
      set(algos "Single1D" "Single2D" "Tubular" "Coupled2D")
      foreach(QUICC_MPIALGO ${algos})
         quicc_add_definition(QUICC_MPIALGO)
      endforeach(QUICC_MPIALGO ${algos})
   endif(NOT QUICC_MPIALGO STREQUAL "Auto")
endif(mpiTest)


###################################################
#---------- MPI COMMUNICATION GROUPING -----------#
###################################################

#
# Choose the type of Serial/MPI parallelisation grouping setup.
# Possible options are: Serial, Auto, Single1D, Single2D, Tubular
#

if(mpiTest)
   set(QUICC_GROUPERS_SERIAL "Equation" "Auto")
   set(QUICC_GROUPERS_SINGLE1D "Equation" "Auto" "Single1D")
   set(QUICC_GROUPERS_SINGLE2D "Equation" "Auto" "Single2D")
   set(QUICC_GROUPERS_TUBULAR "Equation" "Auto" "Single1D" "Single2D" "Transform")
   set(QUICC_GROUPERS_COUPLED2D "Equation" "Auto" "Single1D")
   string(TOUPPER "QUICC_GROUPERS_${QUICC_MPIALGO}" upGrouper)

   quicc_provide_choice(${upGrouper} "Transform grouping" QUICC_TRANSGROUPER groupTest)

   if(groupTest)
      if(QUICC_MPI)
         if(NOT QUICC_TRANSGROUPER STREQUAL "Auto")
            quicc_add_definition(QUICC_TRANSGROUPER)
            if(NOT QUICC_TRANSGROUPER STREQUAL "Equation")
               set(QUICC_TRANSGROUPER "Equation")
               quicc_add_definition(QUICC_TRANSGROUPER)
            endif(NOT QUICC_TRANSGROUPER STREQUAL "Equation")
         else(NOT QUICC_TRANSGROUPER STREQUAL "Auto")
            set(groupers "Equation" "Single1D" "Single2D" "Transform")
            foreach(QUICC_TRANSGROUPER ${groupers})
               quicc_add_definition(QUICC_TRANSGROUPER)
            endforeach(QUICC_TRANSGROUPER ${groupers})
         endif(NOT QUICC_TRANSGROUPER STREQUAL "Auto")
      else(QUICC_MPI)
         set(QUICC_TRANSGROUPER "Equation")
         quicc_add_definition(QUICC_TRANSGROUPER)
      endif(QUICC_MPI)
   endif(groupTest)
endif(mpiTest)


###################################################
#--------------- MPI DATA PACKING ----------------#
###################################################

if(QUICC_MPI)
   set(QUICC_MPIPACKS "MPI" "Manual")
   if(QUICC_MPIPACK STREQUAL "")
      set(QUICC_MPIPACK "MPI")
   endif(QUICC_MPIPACK STREQUAL "")
   mark_as_advanced(FORCE QUICC_MPIPACK)
   quicc_provide_choice(QUICC_MPIPACKS "MPI data packing" QUICC_MPIPACK mpipackTest)

   if(mpipackTest)
      quicc_add_definition(QUICC_MPIPACK)
   endif(mpipackTest)
endif(QUICC_MPI)


###################################################
#--------------- MPI COMMUNICATION ---------------#
###################################################

if(QUICC_MPI)
   set(QUICC_MPICOMMS "SendRecv" "AllToAll")
   if(QUICC_MPICOMM STREQUAL "")
      set(QUICC_MPICOMM "AllToAll")
   endif(QUICC_MPICOMM STREQUAL "")
   mark_as_advanced(FORCE QUICC_MPICOMM)
   quicc_provide_choice(QUICC_MPICOMMS "MPI communication" QUICC_MPICOMM mpicommTest)

   if(mpicommTest)
      quicc_add_definition(QUICC_MPICOMM)
   endif(mpicommTest)
endif(QUICC_MPI)


###################################################
#------------- BOUNDARY CONDITIONS ---------------#
###################################################

#
# Choose the type of boundray method the code is using.
# Possible options are: Galerkin, Tau
#
set(QUICC_BOUNDARYMETHODS "Galerkin" "Tau")

quicc_provide_choice(QUICC_BOUNDARYMETHODS "Boundary method" QUICC_BOUNDARYMETHOD boundaryTest)

if(boundaryTest)
   quicc_add_definition(QUICC_BOUNDARYMETHOD)
endif(boundaryTest)


###################################################
#--------------- TIME INTEGRATORS ----------------#
###################################################

#
# Choose the type of time integration the code is using.
# Possible options are: ImExRKCB2, ImExRKCB3a, ImExRKCB3b, ImExRKCB3c, ImExRKCB3d, ImExRKCB3e, ImExRKCB3f, ImExRKCB4, ImExRK3, ImExSBDF2
#
set(QUICC_TIMESTEPPERS "ImExRKCB2" "ImExRKCB3a" "ImExRKCB3b" "ImExRKCB3c" "ImExRKCB3d" "ImExRKCB3e" "ImExRKCB3f" "ImExRKCB4" "ImExRK3" "ImExSBDF2")

quicc_provide_choice(QUICC_TIMESTEPPERS "Time integrator" QUICC_TIMESTEPPER timeTest)

if(timeTest)
   quicc_add_definition(QUICC_TIMESTEPPER)
endif(timeTest)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%----------------------- DIRECTORY STRUCTURE ---------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


###################################################
#--------- SOURCE DIRECTORIES AND FILES ----------#
###################################################

#
# List of source subdirectories for the Framework part
#
set(MHDFrameworkSrcDirs
   Base
   Exceptions
   Framework
   Timers
)

##### DEVELOPMENT: NEEDS TO BE REMOVED ONCE CODE IS READY
message(STATUS "!!!WARNING!!! OVERRIDING ACTUAL Framework.cmake DIRECTORY LIST WITH NOTHING")
set(MHDFrameworkSrcDirs )

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR} MHDFrameworkSrcDirs)
