# Create list of sources
set(MHDSources
   )

if(QUICC_SPATIALSCHEME STREQUAL "T")
   list(APPEND MHDSources 
      TScheme.cpp
      IRegular1DScheme.cpp
      )
elseif(QUICC_SPATIALSCHEME STREQUAL "F")
   list(APPEND MHDSources 
      FScheme.cpp
      IRegular1DScheme.cpp
      )
endif(QUICC_SPATIALSCHEME STREQUAL "T")
