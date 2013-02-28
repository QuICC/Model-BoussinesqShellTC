/** \file TimerTools.hpp
 *  \brief Implementation of timer tools
 *
 *  \mhdBug Needs test
 */

#ifndef TIMERTOOLS_HPP
#define TIMERTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Timers/ITimer.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of timer tools
    */
   class TimerTools
   {
      public:
         /**
          * @brief Reset timer (stop and restart) and return elapsed time
          */
         static MHDFloat reset(ITimer& rTimer);
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         TimerTools();

         /**
          * @brief Destructor
          */
         ~TimerTools();

   };

}

#endif // TIMERTOOLS_HPP
