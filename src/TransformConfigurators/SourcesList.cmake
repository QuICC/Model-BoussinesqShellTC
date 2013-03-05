set(MHDSources
   BackwardConfigurator.cpp
   ForwardConfigurator.cpp
   BackwardSerialConfigurator.cpp
   ForwardSerialConfigurator.cpp
)

# Add in the MPI based configuration if required
if(GEOMHDISCC_MPI)
   list(APPEND MHDSources
      BackwardSingle1DConfigurator.cpp
      BackwardSingle2DConfigurator.cpp
      BackwardTubularConfigurator.cpp
      ForwardSingle1DConfigurator.cpp
      ForwardSingle2DConfigurator.cpp
      ForwardTubularConfigurator.cpp
   )
endif(GEOMHDISCC_MPI)
