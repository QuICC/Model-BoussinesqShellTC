#
# Include function library
#
include(../cmake.d/functions.cmake)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%-------------------------- CONFIGURABLE -------------------------------%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


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
set(MHDAllSrcDirs
   Base
   Debug
   Communicators
   Diagnostics
   Enums
   Equations
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
   PhysicalModels
   PhysicalOperators
   PolynomialTransforms
   Profiler
   Python
   Quadratures
   Resolutions
   Simulation
   SparseSolvers
   SpatialSchemes
   StorageProfiler
   Timers
   Timesteppers
   TransformCoordinators
   TransformConfigurators
   TransformGroupers
   TypeSelectors
   Variables
   ../External/Interfaces
)


geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR} MHDAllSrcDirs)
