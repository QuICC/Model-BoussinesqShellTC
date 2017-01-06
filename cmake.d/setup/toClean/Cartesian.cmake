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
message(STATUS "************** Cartesian setup ****************")
message(STATUS "***********************************************")


###################################################
#---------------- SPATIAL SCHEME -----------------#
###################################################

#
# Choose the type of spatial scheme.
# Possible options are: 
#
if(GEOMHDISCC_CODEDIM STREQUAL "1D")
   set(GEOMHDISCC_CARTESIAN_SCHEMES "F" "T")
endif(GEOMHDISCC_CODEDIM STREQUAL "1D")
if(GEOMHDISCC_CODEDIM STREQUAL "2D")
   set(GEOMHDISCC_CARTESIAN_SCHEMES "FF" "TF" "FT" "TT")
endif(GEOMHDISCC_CODEDIM STREQUAL "2D")
if(GEOMHDISCC_CODEDIM STREQUAL "3D")
   set(GEOMHDISCC_CARTESIAN_SCHEMES "FFF" "TFF" "FTF" "FFT" "TTF" "TFT" "FTT" "TTT")
endif(GEOMHDISCC_CODEDIM STREQUAL "3D")

quicc_provide_choice(GEOMHDISCC_CARTESIAN_SCHEMES "Spatial scheme" GEOMHDISCC_SPATIALSCHEME schmTest)

option(GEOMHDISCC_FASTCHEBY "Use fast Chebyshev transform?" ON)

if(schmTest)
   if(NOT GEOMHDISCC_FASTCHEBY)
      string(REGEX REPLACE "T" "Tp" tmp ${GEOMHDISCC_SPATIALSCHEME})
      set(GEOMHDISCC_SPATIALSCHEME ${tmp})
      message(STATUS " --> Use slow Chebyshev transform")
   endif(NOT GEOMHDISCC_FASTCHEBY)

   quicc_add_definition(GEOMHDISCC_SPATIALSCHEME)
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
set(MHDCartesianSrcDirs
   Enums
   Equations
   Transforms
)

quicc_append_sources(All_Srcs ${GEOMHDISCC_SRC_DIR}/Cartesian MHDCartesianSrcDirs)
