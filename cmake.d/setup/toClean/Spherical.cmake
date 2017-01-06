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
message(STATUS "************** Spherical setup ****************")
message(STATUS "***********************************************")


###################################################
#---------------- SPATIAL SCHEME -----------------#
###################################################

#
# Choose the type of spatial scheme.
# Possible options are: 
#
if(QUICC_CODEDIM STREQUAL "1D")
   message(ERROR "Can't setup a spherical simulation in 1D!")
endif(QUICC_CODEDIM STREQUAL "1D")
if(QUICC_CODEDIM STREQUAL "2D")
   set(QUICC_SPHERICAL_SCHEMES "SH")
endif(QUICC_CODEDIM STREQUAL "2D")
if(QUICC_CODEDIM STREQUAL "3D")
   set(QUICC_SPHERICAL_SCHEMES "FDSH" "TSH" "WSH")
endif(QUICC_CODEDIM STREQUAL "3D")

quicc_provide_choice(QUICC_SPHERICAL_SCHEMES "Spatial scheme" QUICC_SPATIALSCHEME schmTest)

option(QUICC_FASTCHEBY "Use fast Chebyshev transform?" ON)

if(schmTest)
   if(NOT QUICC_FASTCHEBY)
      string(REGEX REPLACE "T" "Tp" tmp ${QUICC_SPATIALSCHEME})
      set(QUICC_SPATIALSCHEME ${tmp})
      message(STATUS " --> Use slow Chebyshev transform")
   endif(NOT QUICC_FASTCHEBY)

   quicc_add_definition(QUICC_SPATIALSCHEME)
endif(schmTest)


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
set(MHDSphericalSrcDirs
   Enums
   Equations
   Polynomials
   Transforms
)

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR}/Spherical MHDSphericalSrcDirs)
