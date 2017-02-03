# Create list of sources
set(MHDSources
   TimestepCoordinator.cpp
)

# Select Implicit-Explicit RK CB2 scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB2")
   list(APPEND MHDSources 
      ImExRKCB2.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB2")

# Select Implicit-Explicit RK CB3a scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3a")
   list(APPEND MHDSources 
      ImExRKCB3a.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3a")

# Select Implicit-Explicit RK CB3b scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3b")
   list(APPEND MHDSources 
      ImExRKCB3b.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3b")

# Select Implicit-Explicit RK CB3c scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3c")
   list(APPEND MHDSources 
      ImExRKCB3c.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3c")

# Select Implicit-Explicit RK CB3d scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3d")
   list(APPEND MHDSources 
      ImExRKCB3d.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3d")

# Select Implicit-Explicit RK CB3e scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3e")
   list(APPEND MHDSources 
      ImExRKCB3e.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3e")

# Select Implicit-Explicit RK CB3f scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3f")
   list(APPEND MHDSources 
      ImExRKCB3f.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB3f")

# Select Implicit-Explicit RK CB4 scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRKCB4")
   list(APPEND MHDSources 
      ImExRKCB4.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRKCB4")

# Select Implicit-Explicit RK3 scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExRK3")
   list(APPEND MHDSources 
      ImExRK3.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExRK3")

# Select Implicit-Explicit SBDF2 scheme
if(QUICC_TIMESTEPPER STREQUAL "ImExSBDF2")
   list(APPEND MHDSources 
      ImExSBDF2.cpp
      )
endif(QUICC_TIMESTEPPER STREQUAL "ImExSBDF2")
