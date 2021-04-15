
set(MHDSources
)

# Select FFTW transform tests
if(QUICC_FFT STREQUAL "FFTW")
   list(APPEND MHDSources 
      FftwTests.cpp
      )
endif(QUICC_FFT STREQUAL "FFTW")

# Select cuFFT transforms
#if(QUICC_FFT STREQUAL "cuFFT")
#   list(APPEND MHDSources 
#     CuFftTests.cpp
#      )
#endif(QUICC_FFT STREQUAL "cuFFT")
