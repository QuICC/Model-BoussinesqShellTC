if(QUICC_STORAGEPROFILE)
   # Create list of sources
   set(MHDSources
      StorageProfiler.cpp
      StorageProfilerTools.cpp
   )
else(QUICC_STORAGEPROFILE)
   # Create empty list of sources
   set(MHDSources
   )
endif(QUICC_STORAGEPROFILE)
