#
# Project information
#
project(GeoMHDiSCC CXX)

#
# Check for minimal version for CMake
#
cmake_minimum_required(VERSION 2.8.2)

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
message(STATUS "**************** Global setup *****************")
message(STATUS "***********************************************")

###################################################
#--------------------- PLATFORM ------------------#
###################################################

# Load basic setup for platforms
include(cmake.d/platforms/base.cmake)

# Choose platform
geomhdiscc_provide_choice(GEOMHDISCC_PLATFORMS "Platform" GEOMHDISCC_PLATFORM platformTest)

if(platformTest)
   # Load platform
   geomhdiscc_load_platform(GEOMHDISCC_PLATFORM)

   # Create global configurator
   include(cmake.d/setup/global.cmake)

   # Global code structure
   include(cmake.d/structure/global.cmake)

   # External libraries structure
   include(cmake.d/structure/external.cmake)

   # Base setup
   include(cmake.d/setup/Framework.cmake)

   # Debug setup
   include(cmake.d/setup/Debug.cmake)

   # Automatic setup
   include(cmake.d/setup/automatic.cmake)

   # DEBUG STATEMENT TO COMPILE EVERYTHING
   include(cmake.d/setup/AllSrc.cmake)

   # Executables setup
   include(cmake.d/setup/Executables.cmake)

   # TestSuite setup
   include(cmake.d/setup/TestSuite.cmake)

endif(platformTest)
