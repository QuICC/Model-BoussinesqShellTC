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
message(STATUS "**************** Debug setup ******************")
message(STATUS "***********************************************")

###################################################
#------------------- PROFILING -------------------#
###################################################

#
# Use internal profiler and type of profiler
#
option(GEOMHDISCC_PROFILE "Activate internal profiler?" OFF)

if(GEOMHDISCC_PROFILE)
   add_definitions("-DGEOMHDISCC_PROFILE")

   set(GEOMHDISCC_PROFILERS "Coarse" "Detailed")

   geomhdiscc_provide_choice(GEOMHDISCC_PROFILERS "Profiler" GEOMHDISCC_PROFILER profilerTest)

   if(profilerTest)
      geomhdiscc_add_definition(GEOMHDISCC_PROFILER)
   endif(profilerTest)
endif(GEOMHDISCC_PROFILE)

###################################################
#--------------- STORAGE PROFILING ---------------#
###################################################

#
# Used storage requirements profiler?
#
option(GEOMHDISCC_STORAGEPROFILE "Activate internal storage profiler?" OFF)

if(GEOMHDISCC_STORAGEPROFILE)
   add_definitions("-DGEOMHDISCC_STORAGEPROFILE")

   set(GEOMHDISCC_STORAGEPROFILERS "Coarse" "Detailed")

   geomhdiscc_provide_choice(GEOMHDISCC_STORAGEPROFILERS "Storage profiler" GEOMHDISCC_STORAGEPROFILER storageTest)

   if(storageTest)
      geomhdiscc_add_definition(GEOMHDISCC_STORAGEPROFILER)
   endif(storageTest)
endif(GEOMHDISCC_STORAGEPROFILE)


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

if(GEOMHDISCC_PROFILE)
   set(MHDDebugSrcDirs ${MHDDebugSrcDirs} Profiler)
endif(GEOMHDISCC_PROFILE)

if(GEOMHDISCC_STORAGEPROFILE)
   set(MHDDebugSrcDirs ${MHDDebugSrcDirs} Storage)
endif(GEOMHDISCC_STORAGEPROFILE)

geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR}/Debug MHDDebugSrcDirs)
