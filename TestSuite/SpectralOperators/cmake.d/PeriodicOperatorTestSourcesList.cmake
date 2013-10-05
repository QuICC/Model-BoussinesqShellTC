# Create list of sources for test case
set(MHDTestSources
   IoTools/Formatter.cpp
   SpectralOperators/IOperator.cpp
   SpectralOperators/ChebyshevOperator.cpp
   SpectralOperators/PeriodicOperator.cpp
   FastTransforms/FftSetup.cpp
   FastTransforms/FftwTools.cpp
   FastTransforms/FftwLibrary.cpp
   FastTransforms/ChebyshevFftwTransform.cpp
)

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
