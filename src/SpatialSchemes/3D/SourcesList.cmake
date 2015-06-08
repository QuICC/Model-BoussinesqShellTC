# Create list of sources
set(MHDSources
   )

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TTTScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TFTScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TFF")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      TFFScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "FFF")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      FFFScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "CFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      CFTScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "AFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      AFTScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "WFT")
   list(APPEND MHDSources 
      IRegular3DScheme.cpp
      WFTScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL")
   list(APPEND MHDSources 
      IRegularSHlScheme.cpp
      BLFlScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM")
   list(APPEND MHDSources 
      IRegularSHmScheme.cpp
      BLFmScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL")
   list(APPEND MHDSources 
      IRegularSHlScheme.cpp
      SLFlScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources 
      IRegularSHmScheme.cpp
      SLFmScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "WLF")
   list(APPEND MHDSources 
      IRegularSHScheme.cpp
      WLFScheme.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TTT")
