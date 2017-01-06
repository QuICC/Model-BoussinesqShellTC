if(QUICC_DEBUG)
   # Create list of sources
   set(MHDSources
      SerialDebugger.cpp
   )
else(QUICC_DEBUG)
   # Create empty list of sources
   set(MHDSources
   )
endif(QUICC_DEBUG)
