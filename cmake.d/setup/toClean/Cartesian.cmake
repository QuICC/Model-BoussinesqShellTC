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
if(QUICC_CODEDIM STREQUAL "1D")
   set(QUICC_CARTESIAN_SCHEMES "F" "T")
endif(QUICC_CODEDIM STREQUAL "1D")
if(QUICC_CODEDIM STREQUAL "2D")
   set(QUICC_CARTESIAN_SCHEMES "FF" "TF" "FT" "TT")
endif(QUICC_CODEDIM STREQUAL "2D")
if(QUICC_CODEDIM STREQUAL "3D")
   set(QUICC_CARTESIAN_SCHEMES "FFF" "TFF" "FTF" "FFT" "TTF" "TFT" "FTT" "TTT")
endif(QUICC_CODEDIM STREQUAL "3D")

quicc_provide_choice(QUICC_CARTESIAN_SCHEMES "Spatial scheme" QUICC_SPATIALSCHEME schmTest)

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
set(MHDCartesianSrcDirs
   Enums
   Equations
   Transforms
)

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR}/Cartesian MHDCartesianSrcDirs)
