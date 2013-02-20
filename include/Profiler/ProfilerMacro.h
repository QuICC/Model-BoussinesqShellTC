/** \file ProfilerMacro.h
 *  \brief Preprocessor macros to setup profiling timers if requested in configuration
 */

#ifndef PROFILERMACRO_H
#define PROFILERMACRO_H

#ifdef GEOMHDISCC_PROFILE
   #ifdef GEOMHDISCC_MPI
      // include MPI profiler
      #include "Profiler/MpiProfiler.hpp"

      namespace GeoMHDiSCC {
         namespace Debug {
            /// Typedef for a profiler
            typedef MpiProfiler  ProfilerMacro;
         }
      }
   #else
      // include serial profiler
      #include "Profiler/SerialProfiler.hpp"

      namespace GeoMHDiSCC {
         namespace Debug {
            /// Typedef for a profiler
            typedef SerialProfiler  ProfilerMacro;
         }
      }
   #endif // GEOMHDISCC_MPI

   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()  ProfilerMacro::init()

   /// Define profiler start macro function
   #define ProfilerMacro_start(P)  ProfilerMacro::start(P)

   /// Define profiler stop macro function
   #define ProfilerMacro_stop(P)  ProfilerMacro::stop(P)

   /// Define profiler printInfo macro function
   #define ProfilerMacro_printInfo()  ProfilerTools::printInfo()

   #ifdef GEOMHDISCC_PROFILER_DETAILED
      /// Define detailed profiler start macro function
      #define DetailedProfilerMacro_start(P)  ProfilerMacro::start(P)

      /// Define detailed profiler stop macro function
      #define DetailedProfilerMacro_stop(P)  ProfilerMacro::stop(P)
   #else
      /// Define empty detailed profiler start macro function
      #define DetailedProfilerMacro_start(P)  

      /// Define empty detailed profiler stop macro function
      #define DetailedProfilerMacro_stop(P)  
   #endif // GEOMHDISCC_PROFILER_DETAILED

#else
   /// Define profiler initialisation macro function
   #define ProfilerMacro_init()

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
