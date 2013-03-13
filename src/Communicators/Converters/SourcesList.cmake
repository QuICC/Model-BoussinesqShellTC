set(MHDSources
)

if(GEOMHDISCC_MPI)
   list(APPEND MHDSources
      MpiConverterTools.cpp
   )
endif(GEOMHDISCC_MPI)
