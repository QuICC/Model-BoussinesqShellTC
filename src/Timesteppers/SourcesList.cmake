# Create list of sources
set(MHDSources
   TimestepCoordinator.cpp
)

# Select Implicit-Explicit RK3 scheme
if(GEOMHDISCC_TIMESTEPPER STREQUAL "ImExRK3")
   list(APPEND MHDSources 
      ImExRK3.cpp
      )
endif(GEOMHDISCC_TIMESTEPPER STREQUAL "ImExRK3")

# Select Implicit-Explicit SBDF2 scheme
if(GEOMHDISCC_TIMESTEPPER STREQUAL "ImExSBDF2")
   list(APPEND MHDSources 
      ImExSBDF2.cpp
      )
endif(GEOMHDISCC_TIMESTEPPER STREQUAL "ImExSBDF2")
