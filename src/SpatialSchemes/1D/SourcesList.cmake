# Create list of sources
set(MHDSources
   )

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "T")
   list(APPEND MHDSources 
      TScheme.cpp
      IRegular1DScheme.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "F")
   list(APPEND MHDSources 
      FScheme.cpp
      IRegular1DScheme.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "T")
