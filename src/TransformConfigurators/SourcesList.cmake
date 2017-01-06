set(MHDSources
   BackwardConfigurator2D.cpp
   ForwardConfigurator2D.cpp
   TransformTreeEdge.cpp
   TransformTree.cpp
   TransformTreeTools.cpp
   TransformPathEdge.cpp
   TransformPath.cpp
)

if(QUICC_SPATIALDIMENSION STREQUAL "3D")
   list(APPEND MHDSources
      BackwardConfigurator3D.cpp
      ForwardConfigurator3D.cpp
      )
endif(QUICC_SPATIALDIMENSION STREQUAL "3D")

if(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      Cartesian3DTransformSteps.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      AnnulusTransformSteps.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      CylinderTransformSteps.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      ShellTransformSteps.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphereTransformSteps.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TF" OR QUICC_SPATIALSCHEME STREQUAL "TT")
   list(APPEND MHDSources
      Cartesian2DTransformSteps.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
