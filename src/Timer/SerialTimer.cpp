/** \file SerialTimer.cpp
 *  \brief Source of the serial timer implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timer/SerialTimer.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SerialTimer::SerialTimer(const bool autostart)
   {
      // Check if timer should be started at creation
      if(autostart)
      {
         this->start();
      }
   }

   void SerialTimer::start()
   {
      // Get the starting timespec
      clock_gettime(CLOCK_REALTIME, &this->mStart);
   }

   void SerialTimer::stop()
   {
      // Get the stopping timespec
      clock_gettime(CLOCK_REALTIME, &this->mStop);
   }

   MHDFloat SerialTimer::time() const
   {
      // return elapsed seconds
      return this->elapsedSeconds();
   }

   MHDFloat SerialTimer::elapsedSeconds() const
   {
      // Compute elapsed seconds between the two timespecs
      return static_cast<MHDFloat>(this->mStop.tv_sec - this->mStart.tv_sec) + static_cast<MHDFloat>(this->mStop.tv_nsec - this->mStart.tv_nsec)/1.0e9;
   }

}
