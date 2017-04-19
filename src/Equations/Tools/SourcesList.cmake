# Create list of sources
set(MHDSources
   IEigenTools.cpp
   EquationConditions.cpp
   EquationSorters.cpp
   EquationTools.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT")
   list(APPEND MHDSources
      NoEigenTools.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFT" OR QUICC_SPATIALSCHEME STREQUAL "AFT" OR QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      Eigen1DTools.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFF")
   list(APPEND MHDSources
      Eigen2DTools.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      Eigen3DTools.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFL")
   list(APPEND MHDSources
      EigenSHlmTools.cpp
      EigenSHlTools.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFM" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      EigenSHlmTools.cpp
      EigenSHmTools.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT")
