# List of valid backends
set(MHDTestBackendFlag  GEOMHDISCC_SPLINALG)

# Create list of sources for test case
set(MHDTestSources
   Base/MathConstants.cpp
   IoTools/Formatter.cpp
   Resolutions/TransformSetup.cpp
   FastTransforms/FftSetup.cpp
   SpectralOperators/IOperator.cpp
   SpectralOperators/ChebyshevOperator.cpp
   SpectralOperators/IBoundary.cpp
   SpectralOperators/ChebyshevBoundary.cpp
   SpectralOperators/BoundaryConditions.cpp
   SpectralOperators/PeriodicOperator.cpp
)

# Add FFTW backend files
if(GEOMHDISCC_FFT STREQUAL "FFTW")
   list(APPEND MHDTestSources 
      FastTransforms/FftwTools.cpp
      FastTransforms/FftwLibrary.cpp
      FastTransforms/ChebyshevFftwTransform.cpp
   )
endif(GEOMHDISCC_FFT STREQUAL "FFTW")

# Include solver specific files
if(GEOMHDISCC_SPLINALG STREQUAL "Pardiso")
   list(APPEND MHDTestSources 
      ../External/Interfaces/Pardiso_Real.cpp
      ../External/Interfaces/Pardiso_Complex.cpp
   )
endif(GEOMHDISCC_SPLINALG STREQUAL "Pardiso")

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
