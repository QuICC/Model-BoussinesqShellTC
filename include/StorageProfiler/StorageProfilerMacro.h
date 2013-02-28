/** \file StorageProfilerMacro.h
 *  \brief Preprocessor macros to setup storage profiling if requested in configuration
 *
 *  \mhdBug Needs test
 */

#ifndef STORAGEPROFILERMACRO_H
#define STORAGEPROFILERMACRO_H

#ifdef GEOMHDISCC_STORAGEPROFILE
   // include storage profiler
   #include "StorageProfiler/MemorySize.hpp"
   #include "StorageProfiler/StorageProfiler.hpp"
   #include "StorageProfiler/StorageProfilerTools.hpp"

   namespace GeoMHDiSCC {
      /// Typedef for a profiler
      typedef Debug::StorageProfiler  StorageProfilerMacro;
   }

   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  Debug::StorageProfiler::update(P,m)

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  Debug::StorageProfilerTools::printInfo()

#else
   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  

#endif // GEOMHDISCC_STORAGEPROFILE

#endif // STORAGEPROFILERMACRO_H
