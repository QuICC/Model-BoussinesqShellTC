# Create list of sources
set(MHDSources
   )

if(QUICC_SPATIALSCHEME STREQUAL "TTT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TTTScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TFTScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TFF")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TFFScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      FFFScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "CFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      CFTScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      AFTScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      WFTScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFL")
   list(APPEND MHDSources 
      IRegularSHlScheme.cpp
      BLFlScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "BLFM")
   list(APPEND MHDSources 
      IRegularSHmScheme.cpp
      BLFmScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFL")
   list(APPEND MHDSources 
      IRegularSHlScheme.cpp
      SLFlScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources 
      IRegularSHmScheme.cpp
      SLFmScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "WLFL")
   list(APPEND MHDSources 
      IRegularSHlScheme.cpp
      WLFlScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "WLFM")
   list(APPEND MHDSources 
      IRegularSHmScheme.cpp
      WLFmScheme.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "TTT")
