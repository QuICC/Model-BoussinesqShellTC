set(MHDSources
   FrameworkBase.cpp
   MpiFramework.cpp     # ONLY DEBUG, SHOULD BE REMOVED
   SerialFramework.cpp  # ONLY DEBUG, SHOULD BE REMOVED
)

# Add in the parallel or serial framework depending on setup
if(GEOMHDISCC_MPI)
   list(APPEND MHDSources MpiFramework.cpp)
else(GEOMHDISCC_MPI)
   list(APPEND MHDSources SerialFramework.cpp)
endif(GEOMHDISCC_MPI)
