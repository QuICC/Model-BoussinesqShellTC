/** 
 * @file TimerTools.cpp
 * @brief Source of the timer tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

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
