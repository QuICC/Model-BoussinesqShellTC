set(MHDSources
   NoIndexConv.cpp
   PMIndexConv.cpp
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL")
   list(APPEND MHDSources
      SHlIndexConv.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM")
   list(APPEND MHDSources
      SHmIndexConv.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL")

if(GEOMHDISCC_MPI)
   list(APPEND MHDSources
      MpiConverterTools.cpp
   )
endif(GEOMHDISCC_MPI)
