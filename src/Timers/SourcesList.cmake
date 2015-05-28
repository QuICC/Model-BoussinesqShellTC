# Create list of sources
set(MHDSources
   ITimer.cpp
   TimerTools.cpp
   ExecutionTimer.cpp
   StageTimer.cpp
)

# Add in the parallel or serial timers depending on setup
if(GEOMHDISCC_MPI)
   list(APPEND MHDSources MpiTimer.cpp)
else(GEOMHDISCC_MPI)
   list(APPEND MHDSources SerialTimer.cpp)
endif(GEOMHDISCC_MPI)
