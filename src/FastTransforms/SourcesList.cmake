set(MHDSources
   FftSetup.cpp
)

# Select FFTW transforms
if(GEOMHDISCC_FFT STREQUAL "FFTW")
   list(APPEND MHDSources 
      FftwLibrary.cpp
      FftwTools.cpp
      FftwTransform.cpp
      ChebyshevFftwTransform.cpp
      SphereChebyshevFftwTransform.cpp
      ShellChebyshevFftwTransform.cpp
      CylinderChebyshevFftwTransform.cpp
      AnnulusChebyshevFftwTransform.cpp
      )
endif(GEOMHDISCC_FFT STREQUAL "FFTW")

# Select cuFFT transforms
if(GEOMHDISCC_FFT STREQUAL "cuFFT")
   list(APPEND MHDSources 
      CuFftLibrary.cpp
      CuFftTools.cpp
      CuFftTransform.cpp
      ChebyshevCuFftTransform.cpp
      )
endif(GEOMHDISCC_FFT STREQUAL "cuFFT")
