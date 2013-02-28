/** \file TimerMacro.h
 *  \brief Preprocessor macros to setup timers depending on parallelisation options
 *
 *  \mhdBug Needs test
 */

#ifndef TIMERMACRO_H
#define TIMERMACRO_H

#ifdef GEOMHDISCC_MPI
   // include MPI timer
   #include "Timers/MpiTimer.hpp"

   namespace GeoMHDiSCC {
      /// Typedef for a generic Timer based on the MpiTimer
      typedef MpiTimer  TimerMacro;
   }
#else
   // include serial timer
   #include "Timers/SerialTimer.hpp"

   namespace GeoMHDiSCC {
      /// Typedef for a generic Timer based on the SerialTimer
      typedef SerialTimer  TimerMacro;
   }
#endif // GEOMHDISCC_MPI

#endif // TIMERMACRO_H
