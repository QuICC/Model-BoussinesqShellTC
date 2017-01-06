# Create list of sources
set(MHDSources
)

if(QUICC_SPATIALSCHEME STREQUAL "FF")
   list(APPEND MHDSources 
      FFScheme.cpp
      IRegular2DScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TF")
   list(APPEND MHDSources 
      TFScheme.cpp
      IRegular2DScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "TT")
   list(APPEND MHDSources 
      TTScheme.cpp
      IRegular2DScheme.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "FF")
