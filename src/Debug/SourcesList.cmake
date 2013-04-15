if(GEOMHDISCC_DEBUG)
   # Create list of sources
   set(MHDSources
      SerialDebugger.cpp
   )
else(GEOMHDISCC_DEBUG)
   # Create empty list of sources
   set(MHDSources
   )
endif(GEOMHDISCC_DEBUG)
