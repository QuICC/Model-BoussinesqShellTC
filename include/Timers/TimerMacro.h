/**
 * @file TimerMacro.h
 * @brief Preprocessor macros to setup timers depending on parallelisation options 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMERMACRO_H
#define TIMERMACRO_H

#ifdef QUICC_MPI
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
#endif // QUICC_MPI

#endif // TIMERMACRO_H
