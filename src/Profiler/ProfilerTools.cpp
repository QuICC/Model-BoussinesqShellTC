/** 
 * @file ProfilerTools.cpp
 * @brief Source of the implementation of tools for the profiling timer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "Base/Typedefs.hpp"
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
         int digits = 2;

         // Create nice looking ouput header
         IoTools::Formatter::printNewline(std::cout);
         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printCentered(std::cout, "Profiling information", '%');
         IoTools::Formatter::printLine(std::cout, '%');

         // Output the timinigs
         std::stringstream oss;
         int base = 27;
         for(int i = 0; i < static_cast<int>(ProfilerMacro::NBREAKPOINT); i++)
         {
            oss << ProfilerMacro::pointName(static_cast<ProfilerMacro::BreakPoint>(i)) << ": " << std::fixed << std::setprecision(digits) << ts(i);

            if(max.size() != 0)
            {
               oss << " / " << std::fixed << std::setprecision(digits) << max(i) << " / " << std::fixed << std::setprecision(digits) << min(i);
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
