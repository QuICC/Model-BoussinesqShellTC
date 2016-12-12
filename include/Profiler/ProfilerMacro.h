/**
 * @file ProfilerMacro.h
 * @brief Preprocessor macros to setup profiling timers if requested in configuration 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROFILERMACRO_H
#define PROFILERMACRO_H

#ifdef GEOMHDISCC_PROFILE
   #include "Profiler/ProfilerTools.hpp"
   #ifdef GEOMHDISCC_MPI
      // include MPI profiler
      #include "Profiler/MpiProfiler.hpp"

      namespace GeoMHDiSCC {
         /// Typedef for a profiler
         typedef Debug::MpiProfiler ProfilerMacro;
      }
   #else
      // include serial profiler
      #include "Profiler/SerialProfiler.hpp"

      namespace GeoMHDiSCC {
         /// Typedef for a profiler
         typedef Debug::SerialProfiler ProfilerMacro;
      }
   #endif // GEOMHDISCC_MPI

   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()  GeoMHDiSCC::ProfilerMacro::init()

   /// Reset profiler reset macro function
   #define ProfilerMacro_reset()  GeoMHDiSCC::ProfilerMacro::reset()

   /// Define profiler start macro function
   #define ProfilerMacro_start(P)  GeoMHDiSCC::ProfilerMacro::start(P)

   /// Define profiler stop macro function
   #define ProfilerMacro_stop(P)  GeoMHDiSCC::ProfilerMacro::stop(P)

   /// Define profiler printInfo macro function
   #define ProfilerMacro_printInfo()  GeoMHDiSCC::Debug::ProfilerTools::printInfo()

   #ifdef GEOMHDISCC_PROFILER_DETAILED
      /// Define detailed profiler start macro function
      #define DetailedProfilerMacro_start(P)  GeoMHDiSCC::ProfilerMacro::start(P)

      /// Define detailed profiler stop macro function
      #define DetailedProfilerMacro_stop(P)  GeoMHDiSCC::ProfilerMacro::stop(P)
   #else
      /// Define empty detailed profiler start macro function
      #define DetailedProfilerMacro_start(P)  

      /// Define empty detailed profiler stop macro function
      #define DetailedProfilerMacro_stop(P)  
   #endif // GEOMHDISCC_PROFILER_DETAILED

#else
   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()

   /// Define profiler reset macro function
   #define ProfilerMacro_reset()

   /// Define empty profiler start macro function
   #define ProfilerMacro_start(P)  

   /// Define empty profiler stop macro function
   #define ProfilerMacro_stop(P)  

   /// Define empty detailed profiler start macro function
   #define DetailedProfilerMacro_start(P)  

   /// Define empty detailed profiler stop macro function
   #define DetailedProfilerMacro_stop(P)  

   /// Define empty profiler printInfo macro function
   #define ProfilerMacro_printInfo()  
#endif // GEOMHDISCC_PROFILE

#endif // PROFILERMACRO_H
