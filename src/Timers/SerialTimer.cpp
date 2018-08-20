/** 
 * @file SerialTimer.cpp
 * @brief Source of the serial timer implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timers/SerialTimer.hpp"

// Project includes
//

namespace QuICC {

   SerialTimer::SerialTimer(const bool autostart)
   {
      // Check if timer should be started at creation
      if(autostart)
      {
         this->start();
      } else
      {
         // Zero initialise start timespec structures
         this->mStart.tv_sec = 0.0;
         this->mStart.tv_nsec = 0.0;

         // Zero initialise stop timespec structures
         this->mStop.tv_sec = 0.0;
         this->mStop.tv_nsec = 0.0;
      }
   }

   SerialTimer::~SerialTimer()
   {
   }

   void SerialTimer::start()
   {
//      // Get the starting timespec
      clock_gettime(CLOCK_REALTIME, &this->mStart);
// Stefano: added on July 2018. Needed for Mac version earlier than 10.12 
//      #ifdef __APPLE__
//      mach_absolute_time();
//      #else
//      clock_gettime(CLOCK_REALTIME, &this->mStart);
//      #endif
   }

   void SerialTimer::stop()
   {
      // Get the stopping timespec
      clock_gettime(CLOCK_REALTIME, &this->mStop);
// Stefano: added on July 2018. Needed for Mac version earlier than 10.12
//      #ifdef __APPLE__
//      mach_absolute_time();
//      #else
//      clock_gettime(CLOCK_REALTIME, &this->mStop);
//      #endif
   }

   MHDFloat SerialTimer::time() const
   {
      // return elapsed seconds
      return this->elapsedSeconds();
   }

   MHDFloat SerialTimer::queryTime() const
   {
      // Store current stop time
      timespec tmp = this->mStop;

      // Get current elasped time
      clock_gettime(CLOCK_REALTIME, &const_cast<SerialTimer *>(this)->mStop);
      MHDFloat current = this->elapsedSeconds();
      const_cast<SerialTimer*>(this)->mStop = tmp;

      // return elapsed seconds
      return current;
   }

   MHDFloat SerialTimer::resetTimer()
   {
      // Stop the timer
      this->stop();

      // Get elapsed time
      MHDFloat tmp = this->time();

      // Set start time to stopping time
      this->mStart = this->mStop;

      // return elapsed time
      return tmp;
   }

   MHDFloat SerialTimer::elapsedSeconds() const
   {
      // Compute elapsed seconds between the two timespecs
      return static_cast<MHDFloat>(this->mStop.tv_sec - this->mStart.tv_sec) + static_cast<MHDFloat>(this->mStop.tv_nsec - this->mStart.tv_nsec)/1.0e9;
   }

}
