###################################################
#----------------- PYTHON SCRIPTS ----------------#
###################################################

#
# Install Python scripts
#
INSTALL(FILES 
   Python/accuracy_cartesian_1d.py
   Python/accuracy_cartesian_2d.py
   Python/accuracy_cartesian_3d.py
   Python/annulus.py
   Python/cartesian_0d.py 
   Python/cartesian_1d.py 
   Python/cartesian_2d.py 
   Python/cartesian_3d.py 
   Python/cartesian_boundary_1d.py 
   Python/cartesian_boundary_2d.py 
   Python/chebyshev_tools.py
   Python/cylinder.py 
   Python/generate_recurrence.py 
   Python/shell.py 
   Python/sphere.py 
   Python/symbolic.py 
   Python/utils.py 
   Python/test_aft_model.py
   Python/test_blf_model.py
   Python/test_cft_model.py
   Python/test_fff_model.py
   Python/test_tff_model.py
   Python/test_tft_model.py
   Python/test_slf_model.py
   Python/test_ttt_model.py
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
