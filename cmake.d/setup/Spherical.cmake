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
if(GEOMHDISCC_CODEDIM STREQUAL "1D")
   message(ERROR "Can't setup a spherical simulation in 1D!")
endif(GEOMHDISCC_CODEDIM STREQUAL "1D")
if(GEOMHDISCC_CODEDIM STREQUAL "2D")
   set(GEOMHDISCC_SPHERICAL_SCHEMES "SH")
endif(GEOMHDISCC_CODEDIM STREQUAL "2D")
if(GEOMHDISCC_CODEDIM STREQUAL "3D")
   set(GEOMHDISCC_SPHERICAL_SCHEMES "FDSH" "TSH" "WSH")
endif(GEOMHDISCC_CODEDIM STREQUAL "3D")

geomhdiscc_provide_choice(GEOMHDISCC_SPHERICAL_SCHEMES "Spatial scheme" GEOMHDISCC_SPATIALSCHEME schmTest)

option(GEOMHDISCC_FASTCHEBY "Use fast Chebyshev transform?" ON)

if(schmTest)
   if(NOT GEOMHDISCC_FASTCHEBY)
      string(REGEX REPLACE "T" "Tp" tmp ${GEOMHDISCC_SPATIALSCHEME})
      set(GEOMHDISCC_SPATIALSCHEME ${tmp})
      message(STATUS " --> Use slow Chebyshev transform")
   endif(NOT GEOMHDISCC_FASTCHEBY)

   geomhdiscc_add_definition(GEOMHDISCC_SPATIALSCHEME)
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
   SpectralOperators
   Transforms
)

geomhdiscc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR}/Spherical MHDSphericalSrcDirs)
