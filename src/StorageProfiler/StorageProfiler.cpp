/** \file StorageProfiler.cpp
 *  \brief Source of the implementation of a storage profiler
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
#include "StorageProfiler/StorageProfiler.hpp"

// Project includes
//
#include "Base/IO/Ascii/FormatToolbox.hpp"

namespace GeoMHDiSCC {

namespace Debug {

   std::map<StorageProfiler::StoragePoint, MHDFloat> StorageProfiler::requirements = std::map<StorageProfiler::StoragePoint, MHDFloat>();

   void StorageProfiler::update(StorageProfiler::StoragePoint point, MHDFloat memory)
   {
      // Check if entry already exists
      if(StorageProfiler::requirements.count(point))
      {
         // Increment required storage
         StorageProfiler::requirements.at(point) += memory;
      } else
      {
         // Add memory point to required storage
         StorageProfiler::requirements.insert(std::make_pair(point, memory));
      }
   }

   void StorageProfiler::analyze(Array& min, Array& max)
   {
      // Get the "global" requirements from MPI code
      #ifdef GEOMHDISCC_MPI
         // Resize the storage
         Array local = Array::Zero(NMAXSTORAGEPOINT);
         std::map<StoragePoint, MHDFloat>::iterator   reqIt;
         for(reqIt = StorageProfiler::requirements.begin(); reqIt != requirements.end(); reqIt++)
         {
            local(reqIt->first) = reqIt->second;
         }
         min.resize(local.size());
         max.resize(local.size());

         // Get the max values
         MPI_Allreduce(local.data(), max.data(), local.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         // Get the min values
         MPI_Allreduce(local.data(), min.data(), local.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

         // Get the mean values
         MPI_Allreduce(MPI_IN_PLACE, local.data(), local.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         int size;
         // Get the number of CPUs involved
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         // Compute mean requirements per CPU
         local /= static_cast<MHDFloat>(size);
         for(reqIt = StorageProfiler::requirements.begin(); reqIt != StorageProfiler::requirements.end(); reqIt++)
         {
            reqIt->second = local(reqIt->first);
         }
      #endif // GEOMHDISCC_MPI
   }

   void StorageProfiler::createNameMap(std::map<StorageProfiler::StoragePoint, std::string>& map)
   {
      map.insert(std::make_pair(VARIABLES, "Variables"));

      map.insert(std::make_pair(TRANSFORMS, "Transforms"));

      map.insert(std::make_pair(TEMPORARIES, "Temporaries"));

      map.insert(std::make_pair(TIMESTEP, "Timestep"));

      map.insert(std::make_pair(MPI, "MPI"));

      map.insert(std::make_pair(IO, "IO"));

      map.insert(std::make_pair(SCALAR1, "Scalar field 1"));

      map.insert(std::make_pair(SCALAR2, "Scalar field 2"));

      map.insert(std::make_pair(SCALAR3, "Scalar field 3"));

      map.insert(std::make_pair(VECTOR1, "Vector field 1"));

      map.insert(std::make_pair(VECTOR2, "Vector field 2"));

      map.insert(std::make_pair(VECTOR3, "Vector field 3"));

      map.insert(std::make_pair(TRANSFORM1D, "Transform 1D"));

      map.insert(std::make_pair(TRANSFORM2D, "Transform 2D"));

      map.insert(std::make_pair(TRANSFORM3D, "Transform 3D"));

      map.insert(std::make_pair(TRANSFORM1DTEMP, "Transform 1D tmp"));

      map.insert(std::make_pair(TRANSFORM2DTEMP, "Transform 2D tmp"));

      map.insert(std::make_pair(TRANSFORM3DTEMP, "Transform 3D tmp"));

      map.insert(std::make_pair(MPIBUFFERS, "MPI buffers"));

      map.insert(std::make_pair(MPITYPES, "MPI types"));

      map.insert(std::make_pair(MPICOMM, "MPI communication"));
   }

   void StorageProfiler::printInfo()
   {
      // Analyze the data
      Array min;
      Array max;
      StorageProfiler::analyze(min, max);

      std::map<StorageProfiler::StoragePoint, std::string> pointMap;
      StorageProfiler::createNameMap(pointMap);

      // Create nice looking ouput header
      FormatToolbox::printNewline();
      FormatToolbox::printLine('%');
      FormatToolbox::printCentered("Storage profiling information", '%');
      FormatToolbox::printLine('%');

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
      std::map<StoragePoint, MHDFloat>::iterator   reqIt;
      for(reqIt = requirements.begin(); reqIt != requirements.end(); reqIt++)
      {
         // Increment total memory usage
         if(reqIt->first < StorageProfiler::SCALAR1)
         {
            total(0) += reqIt->second;
         }

         // Set unit
         StorageProfiler::setUnit(reqIt->second, memExt, memFactor);

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
         FormatToolbox::printCentered(oss.str(), ' ', base);
         oss.str("");
      }

      // Set unit
      StorageProfiler::setUnit(total(0), memExt, memFactor);

      // Output maximal usage
      oss << "Total" << ": " << std::fixed << std::setprecision(1) << total(0)/memFactor;
      if(max.size() != 0)
      {
         oss << " / " << std::fixed << std::setprecision(1) << total(2)/memFactor << " / " << std::fixed << std::setprecision(1) << total(1)/memFactor;
      }

      // Add size extension
      oss << memExt;
      FormatToolbox::printCentered(oss.str(), ' ', base);
      oss.str("");

      FormatToolbox::printLine('%');
   }

   void StorageProfiler::setUnit(MHDFloat value, std::string &ext, MHDFloat &factor)
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
