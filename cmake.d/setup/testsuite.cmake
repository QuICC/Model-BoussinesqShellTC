###################################################
#----------- TESTSUITE CONFIGURATION -------------#
###################################################

include("cmake.d/testsuite/Base.cmake")

if(GEOMHDISCC_CODEGEOM STREQUAL "Cartesian")
   include("cmake.d/testsuite/Cartesian.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Cartesian")

if(GEOMHDISCC_CODEGEOM STREQUAL "Cylindrical")
   include("cmake.d/testsuite/Cylindrical.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Cylindrical")

if(GEOMHDISCC_CODEGEOM STREQUAL "Spherical")
   include("cmake.d/testsuite/Spherical.cmake")
endif(GEOMHDISCC_CODEGEOM STREQUAL "Spherical")
