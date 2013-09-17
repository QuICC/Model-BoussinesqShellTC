/** 
 * @file StorageProfilerTools.cpp
 * @brief Source of the implementation of tools for the storage profiler
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
#include "StorageProfiler/StorageProfilerTools.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "IoTools/Formatter.hpp"
#include "StorageProfiler/StorageProfiler.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   void StorageProfilerTools::printInfo()
   {
      // Analyze the data
      Array min;
      Array max;
      StorageProfiler::analyze(min, max);

      std::map<StorageProfiler::StoragePoint, std::string> pointMap;
      StorageProfiler::createNameMap(pointMap);

      // Create nice looking ouput header
      IoTools::Formatter::printNewline(std::cout);
      IoTools::Formatter::printLine(std::cout, '%');
      IoTools::Formatter::printCentered(std::cout, "Storage profiling information", '%');
      IoTools::Formatter::printLine(std::cout, '%');

      // Output the storage profile
      std::stringstream oss;
      Array total = Array::Zero(1);
      if(max.size() != 0)
      {
         total = Array::Zero(min.size());
      }
      int base = 27;
      std::string memExt;
      MHDFloat memFactor = 0.0;
      std::map<StorageProfiler::StoragePoint, MHDFloat>::iterator   reqIt;
      for(reqIt = StorageProfiler::requirements.begin(); reqIt != StorageProfiler::requirements.end(); reqIt++)
      {
         // Increment total memory usage
         if(reqIt->first < StorageProfiler::SCALAR1)
         {
            total(0) += reqIt->second;
         }

         // Set unit
         StorageProfilerTools::setUnit(reqIt->second, memExt, memFactor);

         oss << pointMap.at(reqIt->first) << ": " << std::fixed << std::setprecision(1) << reqIt->second/memFactor;
         if(max.size() != 0)
         {
            if(reqIt->first < StorageProfiler::SCALAR1)
            {
               total(1) += min(reqIt->first);
               total(2) += max(reqIt->first);
            }

            oss << " / " << std::fixed << std::setprecision(1) << max(reqIt->first)/memFactor << " / " << std::fixed << std::setprecision(1) << min(reqIt->first)/memFactor;
         }

         // Add size extension
         oss << memExt;
         IoTools::Formatter::printCentered(std::cout, oss.str(), ' ', base);
         oss.str("");
      }

      // Set unit
      StorageProfilerTools::setUnit(total(0), memExt, memFactor);

      // Output maximal usage
      oss << "Total" << ": " << std::fixed << std::setprecision(1) << total(0)/memFactor;
      if(max.size() != 0)
      {
         oss << " / " << std::fixed << std::setprecision(1) << total(2)/memFactor << " / " << std::fixed << std::setprecision(1) << total(1)/memFactor;
      }

      // Add size extension
      oss << memExt;
      IoTools::Formatter::printCentered(std::cout, oss.str(), ' ', base);
      oss.str("");

      IoTools::Formatter::printLine(std::cout, '%');
   }

   void StorageProfilerTools::setUnit(MHDFloat value, std::string &ext, MHDFloat &factor)
   {
      if(value < 1024.)
      {
         ext = " B";
         factor = 1.0;
      } else if(value < 1024.*1024.)
      {
         ext = " k";
         factor = 1024.;
      } else if(value < 1024.*1024.*1024.)
      {
         ext = " M";
         factor = 1024.*1024.;
      } else if(value < 1024.*1024.*1024.*1024.)
      {
         ext = " G";
         factor = 1024.*1024.*1024.;
      }
   }

}
}
