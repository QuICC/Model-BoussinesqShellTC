###################################################
#----------------- PYTHON SCRIPTS ----------------#
###################################################

#
# Install Python package and scripts
#
INSTALL(DIRECTORY Python/quicc DESTINATION ${PROJECT_BINARY_DIR}/Python FILES_MATCHING PATTERN "*.py")
INSTALL(DIRECTORY Python/scripts DESTINATION ${PROJECT_BINARY_DIR}/Python FILES_MATCHING PATTERN "*.py")

# Set path to installed python package and scripts
set(GEOMHDISCC_PYTHON_DIR "${PROJECT_BINARY_DIR}/Python")

#
# Configure the Python embedding wrapper
#
configure_file(
   "${PROJECT_SOURCE_DIR}/include/Python/PythonConfig.hpp.in"
   "${PROJECT_BINARY_DIR}/include/Python/PythonConfig.hpp"
   )

# Add binary directory to include search path
include_directories(${PROJECT_BINARY_DIR}/include)
