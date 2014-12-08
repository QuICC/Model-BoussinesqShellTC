# Create list of sources
set(MHDSources
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "FF")
   list(APPEND MHDSources 
      FFScheme.cpp
      IRegular2DScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TF")
   list(APPEND MHDSources 
      TFScheme.cpp
      IRegular2DScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "TT")
   list(APPEND MHDSources 
      TTScheme.cpp
      IRegular2DScheme.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "FF")
