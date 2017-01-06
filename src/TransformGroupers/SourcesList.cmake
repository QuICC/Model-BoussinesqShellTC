set(MHDSources
   IForwardGrouper2D.cpp
   IBackwardGrouper2D.cpp
)
if(QUICC_SPATIALDIMENSION STREQUAL "3D")
   list(APPEND MHDSources
      IForwardGrouper3D.cpp
      IBackwardGrouper3D.cpp
      )
endif(QUICC_SPATIALDIMENSION STREQUAL "3D")
