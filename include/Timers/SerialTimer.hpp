/**
 * @file SerialTimer.hpp
 * @brief Implementation of a serial timer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SERIALTIMER_HPP
#define SERIALTIMER_HPP

// System includes
//
#include <time.h>

// External includes
//
//#include <time.h>
// Stefano: added in July 2018 for Mac (versions earlier than 10.12)
//#if defined(__APPLE__)
//#include <mach/mach.h>
//#include <mach/mach_time.h>
//#endif
// Project includes
//
#include "Timers/ITimer.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a serial timer
    */
   class SerialTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit SerialTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~SerialTimer();

         /**
          * @brief Start clock
          */
         virtual void start();

         /**
          * @brief Stop clock
          */
         virtual void stop();

         /**
          * @brief Get elapsed time
          */
         virtual MHDFloat time() const;

         /**
          * @brief Query current time without changes to state
          */
         virtual MHDFloat queryTime() const;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual MHDFloat resetTimer();
         
      protected:

      private:
         /**
          * @brief Start
          */
         timespec mStart;

         /**
          * @brief Stop
          */
         timespec mStop;

         /**
          * @brief Compute elapsed seconds between start and stop
          */
         MHDFloat elapsedSeconds() const;
   };

}

#endif // SERIALTIMER_HPP
