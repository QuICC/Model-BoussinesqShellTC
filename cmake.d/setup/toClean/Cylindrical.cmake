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
message(STATUS "************** Cylindrical setup **************")
message(STATUS "***********************************************")


###################################################
#---------------- SPATIAL SCHEME -----------------#
###################################################

#
# Choose the type of spatial scheme.
# Possible options are: 
#
if(QUICC_CODEDIM STREQUAL "1D")
   message(ERROR "Can't setup a cylindrical simulation in 1D!")
endif(QUICC_CODEDIM STREQUAL "1D")
if(QUICC_CODEDIM STREQUAL "2D")
   set(QUICC_CYLINDRICAL_SCHEMES "TF" "WF")
endif(QUICC_CODEDIM STREQUAL "2D")
if(QUICC_CODEDIM STREQUAL "3D")
   set(QUICC_CYLINDRICAL_SCHEMES "TFT" "TFF" "WFT" "WFF")
endif(QUICC_CODEDIM STREQUAL "3D")

quicc_provide_choice(QUICC_CYLINDRICAL_SCHEMES "Spatial scheme" QUICC_SPATIALSCHEME schmTest)

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
set(MHDCylindricalSrcDirs
   Enums
   Equations
   Polynomials
   Transforms
)

quicc_append_sources(All_Srcs ${QUICC_SRC_DIR}/Cylindrical MHDCylindricalSrcDirs)
