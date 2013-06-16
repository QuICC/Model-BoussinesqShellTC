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

#
# List of source subdirectories for the Framework part
#
set(MHDAllSrcDirs
   Base
   Debug
   Communicators
   Diagnostics
   EigenSolver
   Equations
   Enums
   Exceptions
   FastTransforms
   Framework
   Generator
   IoAscii
   IoBinary
   IoConfig
   IoControl
   IoHdf5
   IoTools
   IoVariable
   IoXml
   LoadSplitter
   Profiler
   Resolutions
   Simulation
   SparseSolvers
   SpatialSchemes
   SpectralOperators
   StorageProfiler
   Timers
   Timesteppers
   TransformConfigurators
   TransformGroupers
   TypeSelectors
   Variables
)


geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR} MHDAllSrcDirs)
