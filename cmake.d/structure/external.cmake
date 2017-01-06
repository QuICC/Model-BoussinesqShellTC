###################################################
#-------------- EXTERNAL DIRECTORIES -------------#
###################################################

set(QUICC_EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/External)
set(QUICC_EIGEN_DIR ${QUICC_EXTERNAL_DIR}/eigen3)
set(QUICC_XML_DIR ${QUICC_EXTERNAL_DIR}/rapidxml)
if(QUICC_MULTPRECISION)
   set(QUICC_MPWRAPPER_DIR ${QUICC_EXTERNAL_DIR}/mpfrc++)
endif(QUICC_MULTPRECISION)

###################################################
#-------- EXTERNAL INCLUDE DIRECTORIES -----------#
###################################################

#
# Set External include directories
#
include_directories(${QUICC_EIGEN_DIR})
include_directories(${QUICC_EIGEN_DIR}/unsupported)
include_directories(${QUICC_XML_DIR})
if(QUICC_MULTPRECISION)
   include_directories(${QUICC_MPWRAPPER_DIR})
endif(QUICC_MULTPRECISION)
