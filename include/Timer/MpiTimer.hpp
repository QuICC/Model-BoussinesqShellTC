/** \file MpiTimer.hpp
 *  \brief Implementation of a MPI timer
 */

#ifndef MPITIMER_HPP
#define MPITIMER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Timer/ITimer.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of a MPI timer
    */
   class MpiTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         MpiTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~MpiTimer() {};

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
         virtual MHDFloat resetTimer();
         
      protected:

      private:
         /**
          * @brief Start
          */
         MHDFloat mStart;

         /**
          * @brief Stop
          */
         MHDFloat mStop;
   };

}

#endif // MPITIMER_HPP
