/** \file DebuggerMacro.h
 *  \brief Preprocessor macros to setup debug calls if requested in configuration
 */

#ifndef DEBUGGERMACRO_H
#define DEBUGGERMACRO_H

#ifdef GEOMHDISCC_DEBUG
   #include "Debug/SerialDebugger.hpp"

   namespace GeoMHDiSCC {
      /// Typedef for a profiler
      typedef Debug::SerialDebugger DebuggerMacro;
   }

   /// Define debug enter macro function
   #define DebuggerMacro_enter(P)  DebuggerMacro::enter(P)

   /// Define debug leave macro function
   #define DebuggerMacro_leave(P)  DebuggerMacro::leave(P)

#else
   /// Define empty debug enter  macro function
   #define DebuggerMacro_enter(P)  

   /// Define empty debug leave macro function
   #define DebuggerMacro_leave(P)  
#endif // GEOMHDISCC_DEBUG

#endif // DEBUGGERMACRO_H
