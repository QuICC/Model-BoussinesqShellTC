set(MHDSources
   BackwardConfigurator.cpp
   ForwardConfigurator.cpp
   ForwardSerialConfigurator.cpp
   ProjectorBranch.cpp
   ProjectorTree.cpp
   ProjectorTreeTools.cpp
   IntegratorBranch.cpp
   IntegratorTree.cpp
   IntegratorTreeTools.cpp
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      CartesianTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "CFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      CylindricalTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFL")
   list(APPEND MHDSources
      SphericalTransformSteps.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")

# Add in the MPI based configuration if required
if(GEOMHDISCC_MPI)
   if(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
      list(APPEND MHDSources
         ForwardSingle1DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single1D" OR GEOMHDISCC_MPIALGO STREQUAL "Fixed")
   if(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
      list(APPEND MHDSources
         ForwardSingle2DConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Single2D")
   if(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
      list(APPEND MHDSources
         ForwardTubularConfigurator.cpp
      )
   endif(GEOMHDISCC_MPIALGO STREQUAL "Tubular")
endif(GEOMHDISCC_MPI)
