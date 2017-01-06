if(QUICC_PROFILE)
   # Create list of sources
   set(MHDSources
      ProfilerBase.cpp
      ProfilerTools.cpp
   )

   # Add in the parallel or serial profiler depending on setup
   if(QUICC_MPI)
      list(APPEND MHDSources MpiProfiler.cpp)
   else(QUICC_MPI)
      list(APPEND MHDSources SerialProfiler.cpp)
   endif(QUICC_MPI)
else(QUICC_PROFILE)
   # Create empty list of sources
   set(MHDSources
   )
endif(QUICC_PROFILE)
