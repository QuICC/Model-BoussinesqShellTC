/** \file TimerTools.cpp
 *  \brief Source of the timer tools
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timers/TimerTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   MHDFloat TimerTools::reset(ITimer& rTimer)
   {
      // Stop the timer
      rTimer.stop();

      // Get elapsed time
      MHDFloat tmp = rTimer.time();

      // Restart timer
      rTimer.start();

      // return elapsed time
      return tmp;
   }

}