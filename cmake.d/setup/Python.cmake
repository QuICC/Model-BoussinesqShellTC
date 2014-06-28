###################################################
#----------------- PYTHON SCRIPTS ----------------#
###################################################

#
# Install Python scripts
#
INSTALL(DIRECTORY Python/geomhdiscc DESTINATION ${PROJECT_BINARY_DIR}/Python)

# Set path to installed python scripts
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
