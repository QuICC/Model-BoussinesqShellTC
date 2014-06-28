/**
 * @file ITimer.hpp
 * @brief Implementation of a timer interface 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ITIMER_HPP
#define ITIMER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a timer interface
    */
   class ITimer
   {
      public:
         /**
          * @brief Constructor
          */
         ITimer();

         /**
          * @brief Destructor
          */
         virtual ~ITimer();

         /**
          * @brief Start clock
          */
         virtual void start() = 0;

         /**
          * @brief Stop clock
          */
         virtual void stop() = 0;

         /**
          * @brief Get elapsed time
          */
         virtual MHDFloat time() const = 0;

         /**
          * @brief Get current time without stoping timer
          */
         virtual MHDFloat queryTime() const = 0;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual MHDFloat resetTimer() = 0;
         
      protected:

      private:
   };

}

#endif // ITIMER_HPP
