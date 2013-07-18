set(MHDSources
   BackwardConfigurator.cpp
   ForwardConfigurator.cpp
   BackwardSerialConfigurator.cpp
   ForwardSerialConfigurator.cpp
)

# Add in the MPI based configuration if required
if(GEOMHDISCC_MPI)
   if(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
      list(APPEND MHDSources
         BackwardSingle1DConfigurator.cpp
         ForwardSingle1DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
   if(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
      list(APPEND MHDSources
         BackwardSingle2DConfigurator.cpp
         ForwardSingle2DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
   if(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
      list(APPEND MHDSources
         BackwardTubularConfigurator.cpp
         ForwardTubularConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
endif(GEOMHDISCC_MPI)
