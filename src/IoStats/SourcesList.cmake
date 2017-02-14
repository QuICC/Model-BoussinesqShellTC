# Create list of sources
set(MHDSources
   IStatisticsAsciiEWriter.cpp
   AvgTags.cpp
   KurtTags.cpp
   RmsTags.cpp
   SkewTags.cpp
)

if(QUICC_SPATIALSCHEME STREQUAL "TTT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFF")
   list(APPEND MHDSources
      Cartesian1DScalarAvgWriter.cpp
      Cartesian1DScalarKurtWriter.cpp
      Cartesian1DScalarRMSWriter.cpp
      Cartesian1DScalarSkewWriter.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT" OR QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL" OR QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL" OR QUICC_SPATIALSCHEME STREQUAL "BLFM" OR QUICC_SPATIALSCHEME STREQUAL "WLFL" OR QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT")
