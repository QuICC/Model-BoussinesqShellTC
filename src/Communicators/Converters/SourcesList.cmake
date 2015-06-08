set(MHDSources
   NoIndexConv.cpp
   PMIndexConv.cpp
)

if(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL")
   list(APPEND MHDSources
      SHlIndexConv.cpp
      )
elseif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFM" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFM")
   list(APPEND MHDSources
      SHmIndexConv.cpp
      )
endif(GEOMHDISCC_SPATIALSCHEME STREQUAL "SLFL" OR GEOMHDISCC_SPATIALSCHEME STREQUAL "BLFL")

if(GEOMHDISCC_MPI)
   list(APPEND MHDSources
      MpiConverterTools.cpp
   )
endif(GEOMHDISCC_MPI)
