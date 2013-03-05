if(GEOMHDISCC_STORAGEPROFILE)
   # Create list of sources
   set(MHDSources
      StorageProfiler.cpp
      StorageProfilerTools.cpp
   )
else(GEOMHDISCC_STORAGEPROFILE)
   # Create empty list of sources
   set(MHDSources
   )
endif(GEOMHDISCC_STORAGEPROFILE)
