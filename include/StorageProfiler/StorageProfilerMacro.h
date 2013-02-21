/** \file StorageProfilerMacro.h
 *  \brief Preprocessor macros to setup storage profiling if requested in configuration
 */

#ifndef STORAGEPROFILERMACRO_H
#define STORAGEPROFILERMACRO_H

#ifdef GEOMHDISCC_STORAGEPROFILE
   // include storage profiler
   #include "StorageProfiler/MemorySize.hpp"
   #include "StorageProfiler/StorageProfiler.hpp"

   namespace GeoMHDiSCC {
      namespace Debug {
         /// Typedef for a profiler
         typedef StorageProfiler  StorageProfilerMacro;
      }
   }

   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  StorageProfiler::update(P,m)

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  StorageProfiler::printInfo()

#else
   /// Define storage profiler update macro function
   #define StorageProfilerMacro_update(P,m)  

   /// Define storage profiler printInfo macro function
   #define StorageProfilerMacro_printInfo()  

#endif // GEOMHDISCC_STORAGEPROFILE

#endif // STORAGEPROFILERMACRO_H
