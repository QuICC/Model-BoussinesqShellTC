#
# Project information
#
project(GeoMHDiSCC CXX)
include(cmake.d/MacroOutOfSourceBuild.cmake)
MACRO_ENSURE_OUT_OF_SOURCE_BUILD(
   "${CMAKE_PROJECT_NAME} requires an out of source build. Delete CMakeCache.txt and CMakeFiles/"
)

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

   # Python setup
   include(cmake.d/setup/Python.cmake)

   # Scripts setup
   include(cmake.d/setup/Scripts.cmake)

   # Debug setup
   include(cmake.d/setup/Debug.cmake)

   # Automatic setup
   include(cmake.d/setup/automatic.cmake)

   # Executables setup
   include(cmake.d/setup/Executables.cmake)

   # TestSuite setup
   include(cmake.d/setup/TestSuite.cmake)

endif(platformTest)
