set(MHDSources
   FftSetup.cpp
)

# Select FFTW transforms
if(QUICC_FFT STREQUAL "FFTW")
   list(APPEND MHDSources 
      FftwLibrary.cpp
      FftwTools.cpp
      FftwTransform.cpp
      )
   if(QUICC_SPATIALSCHEME STREQUAL "AFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         AnnulusChebyshevFftwTransform.cpp
         )
   elseif(QUICC_SPATIALSCHEME STREQUAL "CFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         CylinderChebyshevFftwTransform.cpp
         ParityTransformTools.cpp
         )
   elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
      list(APPEND MHDSources 
         ShellChebyshevFftwTransform.cpp
         )
   elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM")
      list(APPEND MHDSources 
         SphereChebyshevFftwTransform.cpp
         ParityTransformTools.cpp
         )
   else(QUICC_SPATIALSCHEME STREQUAL "AFT")
      list(APPEND MHDSources 
         ChebyshevFftwTransform.cpp
         )
   endif(QUICC_SPATIALSCHEME STREQUAL "AFT")
endif(QUICC_FFT STREQUAL "FFTW")

# Select cuFFT transforms
if(QUICC_FFT STREQUAL "cuFFT")
   list(APPEND MHDSources 
      CuFftLibrary.cpp
      CuFftTools.cpp
      CuFftTransform.cpp
      ChebyshevCuFftTransform.cpp
      )
endif(QUICC_FFT STREQUAL "cuFFT")
