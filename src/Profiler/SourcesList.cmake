if(GEOMHDISCC_PROFILE)
   # Create list of sources
   set(MHDSources
      ProfilerBase.cpp
      ProfilerTools.cpp
   )

   # Add in the parallel or serial profiler depending on setup
   if(GEOMHDISCC_MPI)
      list(APPEND MHDSources MpiProfiler.cpp)
   else(GEOMHDISCC_MPI)
      list(APPEND MHDSources SerialProfiler.cpp)
   endif(GEOMHDISCC_MPI)
else(GEOMHDISCC_PROFILE)
   # Create empty list of sources
   set(MHDSources
   )
endif(GEOMHDISCC_PROFILE)
