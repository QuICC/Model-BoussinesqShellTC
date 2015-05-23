/** 
 * @file StageTimer.cpp
 * @brief Source of the stage timer implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "Timers/StageTimer.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

   TimerMacro StageTimer::timer = TimerMacro();

   void StageTimer::msg(const std::string& msg)
   {
      IoTools::Formatter::printCentered(std::cout, "(... " + msg + "...)", ' ');
   }

   void StageTimer::start(const std::string& msg)
   {
      StageTimer::timer.start();
      StageTimer::msg(msg);
   }

   void StageTimer::done()
   {
      StageTimer::timer.stop();

      std::stringstream ss;
      ss << std::ceil(10.0*StageTimer::timer.time())/10.0;

      IoTools::Formatter::printCentered(std::cout, "(... done (" + ss.str() + " s) ...)" , ' ');
      IoTools::Formatter::printNewline(std::cout);
   }

   void StageTimer::newStage(const std::string& msg)
   {
      IoTools::Formatter::printLine(std::cout, '-');
      IoTools::Formatter::printCentered(std::cout, "... " + msg + " ...", '*');
      IoTools::Formatter::printLine(std::cout, '-');
      IoTools::Formatter::printNewline(std::cout);
   }

   void StageTimer::completed(const std::string& msg)
   {
      IoTools::Formatter::printLine(std::cout, '-');
      IoTools::Formatter::printCentered(std::cout, msg, '*');
      IoTools::Formatter::printLine(std::cout, '-');
   }

}
