set(MHDSources
   NoIndexConv.cpp
   SHIndexConv.cpp
   PMIndexConv.cpp
)

if(GEOMHDISCC_MPI)
   list(APPEND MHDSources
      MpiConverterTools.cpp
   )
endif(GEOMHDISCC_MPI)
