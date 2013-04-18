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

   /// Define debug message macro function
   #define DebuggerMacro_msg(M,T)  DebuggerMacro::msg(M,T)

   /// Define debug enter macro function
   #define DebuggerMacro_enter(M,T)  DebuggerMacro::enter(M,T)

   /// Define debug leave macro function
   #define DebuggerMacro_leave(M,T)  DebuggerMacro::leave(M,T)

   /// Define debug showValue macro function
   #define DebuggerMacro_showValue(M,T,V)  DebuggerMacro::showValue(M,T,V)

   /// Define debug timer start macro
   #define DebuggerMacro_start(M,T)  DebuggerMacro::msg(std::string("TIMER START: ")+M,T);DebuggerMacro::timer.start()

   /// Define debug timer stop macro
   #define DebuggerMacro_stop(M,T)  DebuggerMacro::timer.stop();DebuggerMacro::showValue(std::string("TIMER STOP: ")+M,T,DebuggerMacro::timer.time())

#else
   /// Define empty message macro function
   #define DebuggerMacro_msg(M,T)

   /// Define empty debug enter  macro function
   #define DebuggerMacro_enter(M,T)  

   /// Define empty debug leave macro function
   #define DebuggerMacro_leave(M,T)  

   /// Define empty showValue macro function
   #define DebuggerMacro_showValue(M,T,V)

   /// Define empty timer start macro
   #define DebuggerMacro_start(M,T)

   /// Define empty timer stop macro
   #define DebuggerMacro_stop(M,T)
#endif // GEOMHDISCC_DEBUG

#endif // DEBUGGERMACRO_H