set(MHDSources
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources 
      PolySetup.cpp
      PolynomialTools.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF")
   list(APPEND MHDSources 
      WorlandTransform.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF")

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF")
   list(APPEND MHDSources 
      AssociatedLegendreTransform.cpp
      AssociatedLegendreFlyTransform.cpp
      AssociatedLegendrePolynomial.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR  GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF")
