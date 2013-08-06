# Create list of sources
set(MHDSources
)

if(GEOMHDISCC_SPLINALG STREQUAL "Pardiso")
   list(APPEND MHDSources Pardiso_Real.cpp Pardiso_Complex.cpp)
endif(GEOMHDISCC_SPLINALG STREQUAL "Pardiso")
