set(MHDSources
   RandomScalarState.cpp
   RandomVectorState.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      CartesianExactScalarState.cpp
      CartesianExactStateIds.cpp
      CartesianExactVectorState.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      AnnulusExactScalarState.cpp
      AnnulusExactStateIds.cpp
      AnnulusExactVectorState.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      CylinderExactScalarState.cpp
      CylinderExactStateIds.cpp
      CylinderExactVectorState.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      ShellExactScalarState.cpp
      ShellExactStateIds.cpp
      ShellExactVectorState.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      SphereExactScalarState.cpp
      SphereExactStateIds.cpp
      SphereExactVectorState.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT" OR QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "TFF" OR QUICC_SPATIALSCHEME STREQUAL "FFF")
