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
#include "IoAscii/DirectAsciiWriter.hpp"

//#define GEOMHDISCC_PERCORE_PROFILING

namespace GeoMHDiSCC {

namespace Debug {

   void ProfilerTools::writeTimings()
   {
      #ifdef GEOMHDISCC_PERCORE_PROFILING
      Array ts;
      ProfilerMacro::getTimings(ts);

      // Initialize output file
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::stringstream ext;
      ext << "." << rank;
      IoAscii::DirectAsciiWriter prof_out("prof_out", ext.str(), "Detailed per core profiling output", "ASCII", "1.0");
      prof_out.initDebug();

      int digits = 3;

      // Create nice looking ouput header
      IoTools::Formatter::printNewline(prof_out.file());
      IoTools::Formatter::printLine(prof_out.file(), '%');
      IoTools::Formatter::printCentered(prof_out.file(), "Profiling information", '%');
      IoTools::Formatter::printLine(prof_out.file(), '%');

      // Output the timinigs
      std::stringstream oss;
      int base = 27;
      for(int i = 0; i < static_cast<int>(ProfilerMacro::NBREAKPOINT); i++)
      {
         // Only print timing if actually used
         if(ts(i) > 0.0)
         {
            oss << ProfilerMacro::pointName(static_cast<ProfilerMacro::BreakPoint>(i)) << ": " << std::fixed << std::setprecision(digits) << ts(i) << " s";
            IoTools::Formatter::printCentered(prof_out.file(), oss.str(), ' ', base);
            oss.str("");
         }
      }

      IoTools::Formatter::printLine(prof_out.file(), '%');
      IoTools::Formatter::printNewline(prof_out.file());

      prof_out.finalizeDebug();
      #endif //GEOMHDISCC_PERCORE_PROFILING
   }

   void ProfilerTools::printInfo()
   {
      #ifdef GEOMHDISCC_PERCORE_PROFILING
      ProfilerTools::writeTimings();
      #endif //GEOMHDISCC_PERCORE_PROFILING

      // Analyze the data
      Array ts;
      Array min;
      Array max;

      ProfilerMacro::analyze(ts, min, max);

      if(FrameworkMacro::allowsIO())
      {
         int digits = 3;

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
            // Only print timing if actually used
            if(ts(i) > 0.0)
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
         }

         IoTools::Formatter::printLine(std::cout, '%');
         IoTools::Formatter::printNewline(std::cout);
      }
   }

}
}
