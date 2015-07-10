set(MHDSources
   BackwardConfigurator2D.cpp
   ForwardConfigurator2D.cpp
   IntegratorBranch2D.cpp
   ProjectorBranch2D.cpp
   ProjectorTree2D.cpp
   IntegratorTree2D.cpp
)

if(GEOMHDISCC_SPATIALDIMENSION STREQUAL "3D")
   list(APPEND MHDSources
      BackwardConfigurator3D.cpp
      ForwardConfigurator3D.cpp
      IntegratorBranch3D.cpp
      ProjectorBranch3D.cpp
      ProjectorTree3D.cpp
      ProjectorTree3DTools.cpp
      IntegratorTree3D.cpp
      IntegratorTree3DTools.cpp
      )
elseif(GEOMHDISCC_SPATIALDIMENSION STREQUAL "2D")
   list(APPEND MHDSources
      ProjectorTree2DTools.cpp
      IntegratorTree2DTools.cpp
      )
endif(GEOMHDISCC_SPATIALDIMENSION STREQUAL "3D")

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      CartesianTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      AnnulusTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "CFT" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      CylinderTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      ShellTransformSteps.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "WFL")
   list(APPEND MHDSources
      SphereTransformSteps.cpp
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
