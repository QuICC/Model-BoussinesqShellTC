set(MHDSources
   IForwardGrouper2D.cpp
   IBackwardGrouper2D.cpp
)
if(GEOMHDISCC_SPATIALDIMENSION STREQUAL "3D")
   list(APPEND MHDSources
      IForwardGrouper3D.cpp
      IBackwardGrouper3D.cpp
      )
endif(GEOMHDISCC_SPATIALDIMENSION STREQUAL "3D")
