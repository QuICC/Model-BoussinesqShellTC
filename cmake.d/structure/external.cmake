###################################################
#-------------- EXTERNAL DIRECTORIES -------------#
###################################################

set(QUICC_EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/External)
set(QUICC_EIGEN_DIR ${QUICC_EXTERNAL_DIR}/eigen3)
set(QUICC_XML_DIR ${QUICC_EXTERNAL_DIR}/rapidxml)

###################################################
#-------- EXTERNAL INCLUDE DIRECTORIES -----------#
###################################################

#
# Set External include directories
#
include_directories(${QUICC_EIGEN_DIR})
include_directories(${QUICC_EIGEN_DIR}/unsupported)
include_directories(${QUICC_XML_DIR})
