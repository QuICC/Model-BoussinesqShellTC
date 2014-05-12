###################################################
#----------------- PYTHON SCRIPTS ----------------#
###################################################

#
# Install Python scripts
#
INSTALL(FILES 
   Python/cartesian.py 
   Python/sphere.py 
   Python/shell.py 
   Python/annulus.py
   Python/cylinder.py 
   Python/utils.py 
   Python/test_tft_model.py
   DESTINATION ${PROJECT_BINARY_DIR}/Python)

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
