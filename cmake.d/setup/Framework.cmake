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
# Possible options are: BigMem, LowMem
#
set(GEOMHDISCC_MEMORYUSAGES "High" "Limited")

geomhdiscc_provide_choice(GEOMHDISCC_MEMORYUSAGES "Memory usage" GEOMHDISCC_MEMORYUSAGE memoryTest)

if(memoryTest)
   geomhdiscc_add_definition(GEOMHDISCC_MEMORYUSAGE)
endif(memoryTest)


###################################################
#------------- MPI PARALLELISATION ---------------#
###################################################

#
# Choose the type of MPI parallelisation or serial setup.
# Possible options are: Serial, Auto, Single1D, Single2D, Tubular
#
set(GEOMHDISCC_MPIALGOS "Serial" "Auto" "Single1D" "Single2D" "Tubular" "Fixed")

geomhdiscc_provide_choice(GEOMHDISCC_MPIALGOS "MPI algorithm" GEOMHDISCC_MPIALGO mpiTest)

if(mpiTest)
   if(NOT GEOMHDISCC_MPIALGO STREQUAL "Serial")
      set(GEOMHDISCC_MPI ON)
      add_definitions("-DGEOMHDISCC_MPI")
   endif(NOT GEOMHDISCC_MPIALGO STREQUAL "Serial")

   if(NOT GEOMHDISCC_MPIALGO STREQUAL "Auto")
      geomhdiscc_add_definition(GEOMHDISCC_MPIALGO)
   else(NOT GEOMHDISCC_MPIALGO STREQUAL "Auto")
      set(algos "Single1D" "Single2D" "Tubular" "Fixed")
      foreach(GEOMHDISCC_MPIALGO ${algos})
         geomhdiscc_add_definition(GEOMHDISCC_MPIALGO)
      endforeach(GEOMHDISCC_MPIALGO ${algos})
   endif(NOT GEOMHDISCC_MPIALGO STREQUAL "Auto")
endif(mpiTest)

###################################################
#---------- MPI COMMUNICATION GROUPING -----------#
###################################################

#
# Choose the type of Serial/MPI parallelisation grouping setup.
# Possible options are: Serial, Auto, Single1D, Single2D, Tubular
#
set(GEOMHDISCC_GROUPERS "Equation" "Auto" "Single1D" "Single2D" "Transform")

geomhdiscc_provide_choice(GEOMHDISCC_GROUPERS "Transform grouping" GEOMHDISCC_TRANSGROUPER groupTest)

if(groupTest)
   if(GEOMHDISCC_MPI)
      if(NOT GEOMHDISCC_TRANSGROUPER STREQUAL "Auto")
         geomhdiscc_add_definition(GEOMHDISCC_TRANSGROUPER)
         if(NOT GEOMHDISCC_TRANSGROUPER STREQUAL "Equation")
            set(GEOMHDISCC_TRANSGROUPER "Equation")
            geomhdiscc_add_definition(GEOMHDISCC_TRANSGROUPER)
         endif(NOT GEOMHDISCC_TRANSGROUPER STREQUAL "Equation")
      else(NOT GEOMHDISCC_TRANSGROUPER STREQUAL "Auto")
         set(groupers "Equation" "Single1D" "Single2D" "Transform")
         foreach(GEOMHDISCC_TRANSGROUPER ${groupers})
            geomhdiscc_add_definition(GEOMHDISCC_TRANSGROUPER)
         endforeach(GEOMHDISCC_TRANSGROUPER ${groupers})
      endif(NOT GEOMHDISCC_TRANSGROUPER STREQUAL "Auto")
   else(GEOMHDISCC_MPI)
      set(GEOMHDISCC_TRANSGROUPER "Equation")
      geomhdiscc_add_definition(GEOMHDISCC_TRANSGROUPER)
   endif(GEOMHDISCC_MPI)
endif(groupTest)


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

geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR} MHDFrameworkSrcDirs)
