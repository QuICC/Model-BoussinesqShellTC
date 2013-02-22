# Create list of sources
set(MHDSources
   ITimer.cpp
   MpiTimer.cpp      # ONLY DEBUG, SHOULD BE REMOVED
   SerialTimer.cpp   # ONLY DEBUG, SHOULD BE REMOVED
   TimerTools.cpp
)

# Add in the parallel or serial timers depending on setup
if(GEOMHDISCC_MPI)
   list(APPEND MHDSources MpiTimer.cpp)
else(GEOMHDISCC_MPI)
   list(APPEND MHDSources SerialTimer.cpp)
endif(GEOMHDISCC_MPI)
