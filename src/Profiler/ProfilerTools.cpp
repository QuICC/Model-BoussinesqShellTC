/** \file ProfilerTools.cpp
 *  \brief Source of the implementation of tools for the profiling timer
 */

// Configuration includes
//
#include "Profiler/ProfilerMacro.h"
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iostream>

// External includes
//

// Class include
//
#include "Profiler/ProfilerTools.hpp"

// Project includes
//
#include "IoTools/Formatter.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   void ProfilerTools::printInfo()
   {
      // Analyze the data
      Array ts;
      Array min;
      Array max;
      ProfilerMacro::analyze(ts, min, max);

      if(FrameworkMacro::allowsIO())
      {
         // Create nice looking ouput header
         IoTools::Formatter::printNewline(std::cout);
         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printCentered(std::cout, "Profiling information", '%');
         IoTools::Formater::printLine(std::cout, '%');

         // Output the timinigs
         std::stringstream oss;
         int base = 27;
         for(int i = 0; i < ProfilerMacro::NBREAKPOINT; i++)
         {
            oss << ProfilerMacro::pointName(static_cast<ProfilerMacro::BreakPoint>(i)) << ": " << std::fixed << std::setprecision(1) << ts(i);

            if(max.size() != 0)
            {
               oss << " / " << std::fixed << std::setprecision(1) << max(i) << " / " << std::fixed << std::setprecision(1) << min(i);
            }

            oss << " s";
            IoTools::Formatter::printCentered(std::cout, oss.str(), ' ', base);
            oss.str("");
         }

         IoTools::Formatter::printLine(std::cout, '%');
      }
   }

}
}
