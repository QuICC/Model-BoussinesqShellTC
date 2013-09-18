# Set backend flag
set(MHDTestBackendFlag GEOMHDISCC_FFT)

# Create list of sources
set(MHDTestSources
   IoTools/Formatter.cpp
   Base/MathConstants.cpp
   Resolutions/TransformSetup.cpp
   FastTransforms/FftSetup.cpp
)

# Add FFTW backend files
if(GEOMHDISCC_FFT STREQUAL "FFTW")
   list(APPEND MHDTestSources 
      FastTransforms/FftwTools.cpp
      FastTransforms/FftwLibrary.cpp
      FastTransforms/FftwTransform.cpp
   )
endif(GEOMHDISCC_FFT STREQUAL "FFTW")

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
