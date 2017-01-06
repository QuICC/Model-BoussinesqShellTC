set(MHDSources
   FrameworkBase.cpp
)

# Add in the parallel or serial framework depending on setup
if(QUICC_MPI)
   list(APPEND MHDSources MpiFramework.cpp)
else(QUICC_MPI)
   list(APPEND MHDSources SerialFramework.cpp)
endif(QUICC_MPI)
