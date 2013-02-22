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
message(STATUS "******* DEVELOPMENT ALL SOURCES SETUP *********")
message(STATUS "***********************************************")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%----------------------- DIRECTORY STRUCTURE ---------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


###################################################
#--------- SOURCE DIRECTORIES AND FILES ----------#
###################################################

message(STATUS "!!!WARNING!!! AllSrc.cmake OVERRIDING ACTUAL DIRECTORY LIST WITH FULL LIST")
add_definitions("-DGEOMHDISCC_PROFILE")
add_definitions("-DGEOMHDISCC_STORAGEPROFILE")

#
# List of source subdirectories for the Framework part
#
set(MHDAllSrcDirs
   Base
   Exceptions
   FastTransforms
   Framework         # SourcesList.txt IS IN DEBUG MODE
   IoAscii
   IoBinary
   #IoHdf5
   IoTools
   IoXml
   #   LoadSplitter
   Profiler
   #   Resolutions
   #   Simulation
   #   SpatialSchemes
   SpectralOperators
   #   StorageProfiler
   Timers            # SourcesList.txt IS IN DEBUG MODE
)


geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR} MHDAllSrcDirs)
