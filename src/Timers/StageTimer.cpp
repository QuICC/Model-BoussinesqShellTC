/** 
 * @file StageTimer.cpp
 * @brief Source of the stage timer implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

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

   StageTimer::StageTimer()
      : TimerMacro(false), mLevel(0)
   {
   }

   StageTimer::~StageTimer()
   {
   }

   void StageTimer::stage(const std::string& msg)
   {
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printHalfLine(std::cout, '/', '\\');
         IoTools::Formatter::printCentered(std::cout, "... " + msg + " ...", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);
      }
   }

   void StageTimer::completed(const std::string& msg)
   {
      if(FrameworkMacro::allowsIO())
      {
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, msg, '*');
         IoTools::Formatter::printHalfLine(std::cout, '\\', '/');
         IoTools::Formatter::printNewline(std::cout);
      }
   }

   void StageTimer::msg(const std::string& msg, const int space)
   {
      if(FrameworkMacro::allowsIO())
      {
         std::cout << std::setfill(' ') << std::setw(space) << "" << msg << std::endl;
      }
   }

   void StageTimer::start(const std::string& msg, const int level)
   {
      this->mLevel = level;

      if(FrameworkMacro::allowsIO())
      {
         StageTimer::msg( "--> " + msg, 8 + this->mLevel*4);
         TimerMacro::start();
      }
   }

   void StageTimer::done()
   {
      if(FrameworkMacro::allowsIO())
      {
         TimerMacro::stop();

         std::stringstream ss;
         ss << std::ceil(10.0*this->time())/10.0;
         StageTimer::msg("    done (" + ss.str() + " s)", 8 + this->mLevel*4);
         IoTools::Formatter::printNewline(std::cout);
      }
   }

}
