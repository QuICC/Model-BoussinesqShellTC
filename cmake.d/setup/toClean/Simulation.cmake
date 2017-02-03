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
set(QUICC_CODEDIM_DIMENSIONS "1D" "2D" "3D")

quicc_provide_choice(QUICC_CODEDIM_DIMENSIONS "Dimension" QUICC_CODEDIM dimTest)

if(dimTest)
   quicc_add_definition(QUICC_CODEDIM)
endif(dimTest)


###################################################
#-------------------- GEOMETRY -------------------#
###################################################

#
# Choose the type of geometry the simulation is contained in.
# Possible options are: Cartesian, Cylindrical, Spherical
#
set(QUICC_CODEGEOM_GEOMETRIES "Cartesian" "Cylindrical" "Spherical")

quicc_provide_choice(QUICC_CODEGEOM_GEOMETRIES "Geometry" QUICC_CODEGEOM geomTest)

if(geomTest)
   quicc_add_definition(QUICC_CODEGEOM)
endif(geomTest)


###################################################
#------------------- PARAMETERS ------------------#
###################################################

#
# Choose the type of nondimensional parameters
# Possible options are: EQRaEm, PrRaXGAs
#
set(QUICC_PARAMETERS_LIST "EQRaEm" "PrRaXGAs")

quicc_provide_choice(QUICC_PARAMETERS_LIST "Parameters" QUICC_PARAMETERS paramTest)

if(paramTest)
   quicc_add_definition(QUICC_PARAMETERS)
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

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR}/Simulation MHDSimulationSrcDirs)

if(QUICC_CODEGEOM STREQUAL "Cartesian")
   include("cmake.d/setup/Cartesian.cmake")
endif(QUICC_CODEGEOM STREQUAL "Cartesian")

if(QUICC_CODEGEOM STREQUAL "Cylindrical")
   include("cmake.d/setup/Cylindrical.cmake")
endif(QUICC_CODEGEOM STREQUAL "Cylindrical")

if(QUICC_CODEGEOM STREQUAL "Spherical")
   include("cmake.d/setup/Spherical.cmake")
endif(QUICC_CODEGEOM STREQUAL "Spherical")
