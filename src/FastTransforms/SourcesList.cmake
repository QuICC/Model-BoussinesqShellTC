set(MHDSources
   FftSetup.cpp
)

# Select FFTW transforms
if(GEOMHDISCC_FFT STREQUAL "FFTW")
   list(APPEND MHDSources 
      FftwLibrary.cpp
      FftwTools.cpp
      FftwTransform.cpp
      )
   if(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         AnnulusChebyshevFftwTransform.cpp
         )
   elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "CFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         CylinderChebyshevFftwTransform.cpp
         ParityTransformTools.cpp
         )
   elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM")
      list(APPEND MHDSources 
         ShellChebyshevFftwTransform.cpp
         )
   elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM")
      list(APPEND MHDSources 
         SphereChebyshevFftwTransform.cpp
         ParityTransformTools.cpp
         )
   else(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         )
   endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
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
