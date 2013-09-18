# Set test name prefix
if(GEOMHDISCC_MPI)
   set(MHDTestPrefix Mpi)
else(GEOMHDISCC_MPI)
   set(MHDTestPrefix Serial)
endif(GEOMHDISCC_MPI)

# Create list of sources for test case
set(MHDTestSources
   IoTools/Formatter.cpp
   Timers/ITimer.cpp
)

if(GEOMHDISCC_MPI)
   list(APPEND MHDTestSources 
      Timers/MpiTimer.cpp
   )
else(GEOMHDISCC_MPI)
   list(APPEND MHDTestSources 
      Timers/SerialTimer.cpp
   )
endif(GEOMHDISCC_MPI)

# Include all files for the framework
include(../src/Framework/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Framework/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the exceptions
include(../src/Exceptions/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Exceptions/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the storage profiler
include(../src/StorageProfiler/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources StorageProfiler/${src})
endforeach(src ${MHDSources})
set(MHDSources )

# Include all files for the profiler
include(../src/Profiler/SourcesList.cmake)
foreach(src ${MHDSources})
   list(APPEND MHDTestSources Profiler/${src})
endforeach(src ${MHDSources})
set(MHDSources )
