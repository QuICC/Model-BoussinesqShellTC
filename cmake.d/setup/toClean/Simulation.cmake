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
message(STATUS "************* Simulation setup ****************")
message(STATUS "***********************************************")

###################################################
#----------------- DIMENSIONALITY ----------------#
###################################################

#
# Choose the dimensionality of the simulation volume.
# Possible options are: 1D, 2D, 3D
#
set(GEOMHDISCC_CODEDIM_DIMENSIONS "1D" "2D" "3D")

quicc_provide_choice(GEOMHDISCC_CODEDIM_DIMENSIONS "Dimension" GEOMHDISCC_CODEDIM dimTest)

if(dimTest)
   quicc_add_definition(GEOMHDISCC_CODEDIM)
endif(dimTest)


###################################################
#-------------------- GEOMETRY -------------------#
###################################################

#
# Choose the type of geometry the simulation is contained in.
# Possible options are: Cartesian, Cylindrical, Spherical
#
set(GEOMHDISCC_CODEGEOM_GEOMETRIES "Cartesian" "Cylindrical" "Spherical")

quicc_provide_choice(GEOMHDISCC_CODEGEOM_GEOMETRIES "Geometry" GEOMHDISCC_CODEGEOM geomTest)

if(geomTest)
   quicc_add_definition(GEOMHDISCC_CODEGEOM)
endif(geomTest)


###################################################
#------------------- PARAMETERS ------------------#
###################################################

#
# Choose the type of nondimensional parameters
# Possible options are: EQRaEm, PrRaXGAs
#
set(GEOMHDISCC_PARAMETERS_LIST "EQRaEm" "PrRaXGAs")

quicc_provide_choice(GEOMHDISCC_PARAMETERS_LIST "Parameters" GEOMHDISCC_PARAMETERS paramTest)

if(paramTest)
   quicc_add_definition(GEOMHDISCC_PARAMETERS)
endif(paramTest)


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
set(MHDSimulationSrcDirs
   General
   PrepMacros
   Enums
   Controls
   Equations
   IO
   System
   Timesteppers
   Variables
)

quicc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR}/Simulation MHDSimulationSrcDirs)

if(GEOMHDISCC_CODEGEOM STREQUAL "Cartesian")
   include("cmake.d/setup/Cartesian.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Cartesian")

if(GEOMHDISCC_CODEGEOM STREQUAL "Cylindrical")
   include("cmake.d/setup/Cylindrical.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Cylindrical")

if(GEOMHDISCC_CODEGEOM STREQUAL "Spherical")
   include("cmake.d/setup/Spherical.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Spherical")
