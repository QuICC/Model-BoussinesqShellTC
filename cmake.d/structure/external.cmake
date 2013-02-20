###################################################
#-------------- EXTERNAL DIRECTORIES -------------#
###################################################

set(GEOMHDISCC_EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/External)
set(GEOMHDISCC_EIGEN_DIR ${GEOMHDISCC_EXTERNAL_DIR}/eigen3)
set(GEOMHDISCC_XML_DIR ${GEOMHDISCC_EXTERNAL_DIR}/rapidxml)
if(GEOMHDISCC_MULTPRECISION)
   set(GEOMHDISCC_MPWRAPPER_DIR ${GEOMHDISCC_EXTERNAL_DIR}/mpfrc++)
endif(GEOMHDISCC_MULTPRECISION)

###################################################
#-------- EXTERNAL INCLUDE DIRECTORIES -----------#
###################################################

#
# Set External include directories
#
include_directories(${GEOMHDISCC_EIGEN_DIR})
include_directories(${GEOMHDISCC_EIGEN_DIR}/unsupported)
include_directories(${GEOMHDISCC_XML_DIR})
if(GEOMHDISCC_MULTPRECISION)
   include_directories(${GEOMHDISCC_MPWRAPPER_DIR})
endif(GEOMHDISCC_MULTPRECISION)
