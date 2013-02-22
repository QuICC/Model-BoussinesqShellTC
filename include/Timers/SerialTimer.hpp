/** \file SerialTimer.hpp
 *  \brief Implementation of a serial timer
 */

#ifndef SERIALTIMER_HPP
#define SERIALTIMER_HPP

// System includes
//
#include <time.h>

// External includes
//

// Project includes
//
#include "Timers/ITimer.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of a serial timer
    */
   class SerialTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         SerialTimer(const bool autostart = false);

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
          * @brief Reset timer (stop and restart)
          */
         virtual MHDFloat resetTimer() = 0;
         
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
