set(MHDSources
   VectorFieldTrivialVisualizer.cpp
   VectorFieldVisualizer.cpp
   ScalarFieldTrivialVisualizer.cpp
   ScalarFieldVisualizer.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      TiltedScalarFieldVisualizer.cpp
      VelocityStreamVisualizer.cpp
      VorticityStreamVisualizer.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      NonlinearScalarFieldVisualizer.cpp
      NonlinearVectorFieldVisualizer.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      SphericalVerticalFieldVisualizer.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphericalVerticalFieldVisualizer.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
