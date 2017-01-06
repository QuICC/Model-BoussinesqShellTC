# Create list of sources
set(MHDSources
   DiagnosticCoordinator.cpp
   ICflWrapper.cpp
   IVelocityWrapper.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      StreamVerticalWrapper.cpp
      CartesianTorPolWrapper.cpp
      CartesianCflWrapper.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM" OR QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphericalTorPolWrapper.cpp
      SphericalCflWrapper.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TF" OR QUICC_SPATIALSCHEME STREQUAL "TT")
   list(APPEND MHDSources
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
