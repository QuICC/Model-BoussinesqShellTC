#
# Include function library
#
include(cmake.d/functions.cmake)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

###################################################
#-------------------- DEBUGGER -------------------#
###################################################

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
   set(QUICC_DEBUG ON)
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if(QUICC_DEBUG)
   message(STATUS "***********************************************")
   message(STATUS "**************** Debug setup ******************")
   message(STATUS "***********************************************")
   message(STATUS " --> Debug: Active")
endif(QUICC_DEBUG)

if(QUICC_DEBUG)
   add_definitions("-DQUICC_DEBUG")
else(QUICC_DEBUG)
   add_definitions("-DQUICC_NO_DEBUG")
endif(QUICC_DEBUG)

###################################################
#------------------- PROFILING -------------------#
###################################################

#
# Use internal profiler and type of profiler
#
option(QUICC_PROFILE "Activate internal profiler?" OFF)

if(NOT QUICC_DEBUG AND QUICC_PROFILE)
   message(STATUS "***********************************************")
   message(STATUS "**************** Debug setup ******************")
   message(STATUS "***********************************************")
endif(NOT QUICC_DEBUG AND QUICC_PROFILE)

if(QUICC_PROFILE)
   add_definitions("-DQUICC_PROFILE")

   set(QUICC_PROFILERS "Coarse" "Detailed")

   quicc_provide_choice(QUICC_PROFILERS "Profiler" QUICC_PROFILER profilerTest)

   if(profilerTest)
      quicc_add_definition(QUICC_PROFILER)
   endif(profilerTest)

   if(QUICC_MPI)
      option(QUICC_PROFILE_PERCORE "Write per core profiling data?" OFF)
      if(QUICC_PROFILE_PERCORE)
         add_definitions("-DQUICC_PROFILE_PERCORE")
      endif(QUICC_PROFILE_PERCORE)
   endif(QUICC_MPI)
endif(QUICC_PROFILE)

###################################################
#--------------- STORAGE PROFILING ---------------#
###################################################

#
# Used storage requirements profiler?
#
option(QUICC_STORAGEPROFILE "Activate internal storage profiler?" OFF)

if(NOT QUICC_DEBUG AND NOT QUICC_PROFILE AND QUICC_STORAGEPROFILE)
   message(STATUS "***********************************************")
   message(STATUS "**************** Debug setup ******************")
   message(STATUS "***********************************************")
endif(NOT QUICC_DEBUG AND NOT QUICC_PROFILE AND QUICC_STORAGEPROFILE)

if(QUICC_STORAGEPROFILE)
   add_definitions("-DQUICC_STORAGEPROFILE")

   set(QUICC_STORAGEPROFILERS "Coarse" "Detailed")

   quicc_provide_choice(QUICC_STORAGEPROFILERS "Storage profiler" QUICC_STORAGEPROFILER storageTest)

   if(storageTest)
      quicc_add_definition(QUICC_STORAGEPROFILER)
   endif(storageTest)
endif(QUICC_STORAGEPROFILE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%----------------------- DIRECTORY STRUCTURE ---------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


###################################################
#--------- SOURCE DIRECTORIES AND FILES ----------#
###################################################

#
# List of source subdirectories for the Base subframework
#

if(QUICC_DEBUG)
   set(MHDDebugSrcDirs ${MHDDebugSrcDirs} Debug)
endif(QUICC_DEBUG)

if(QUICC_PROFILE)
   set(MHDDebugSrcDirs ${MHDDebugSrcDirs} Profiler)
endif(QUICC_PROFILE)

if(QUICC_STORAGEPROFILE)
   set(MHDDebugSrcDirs ${MHDDebugSrcDirs} StorageProfiler)
endif(QUICC_STORAGEPROFILE)

##### DEVELOPMENT: NEEDS TO BE REMOVED ONCE CODE IS READY
message(STATUS "!!!WARNING!!! OVERRIDING ACTUAL Debug.cmake DIRECTORY LIST WITH NOTHING")
set(MHDDebugSrcDirs )

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR} MHDDebugSrcDirs)
