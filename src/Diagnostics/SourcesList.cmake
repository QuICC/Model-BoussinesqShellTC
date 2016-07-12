# Create list of sources
set(MHDSources
   DiagnosticCoordinator.cpp
   ICflWrapper.cpp
   IVelocityWrapper.cpp
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      StreamVerticalWrapper.cpp
      CartesianCflWrapper.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "CFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphericalTorPolWrapper.cpp
      SphericalCflWrapper.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AF")
   list(APPEND MHDSources
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "CF")
   list(APPEND MHDSources
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TT")
   list(APPEND MHDSources
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")
