/** 
 * @file SimulationBoundary.cpp
 * @brief Implementation of a general simulation control structure
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Simulation/SimulationBoundary.hpp"

// Project includes
//
#include "IoTools/HumanToId.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

   SimulationBoundary::SimulationBoundary(const std::map<std::string,int>& bcIds)
   {
      this->convert(bcIds);
   }

   SimulationBoundary::~SimulationBoundary()
   {
   }

   void SimulationBoundary::convert(const std::map<std::string,int>& bcIds)
   {
      std::map<std::string,int>::const_iterator mapIt;

      for(mapIt = bcIds.begin(); mapIt != bcIds.end(); ++mapIt)
      {
         this->mBcs.insert(std::make_pair(IoTools::HumanToId::toPhys(mapIt->first),mapIt->second));
      }

   }

   std::map<std::string,int>  SimulationBoundary::getTagMap() const
   {
      std::map<std::string,int>  tagMap;

      std::map<PhysicalNames::Id,int>::const_iterator mapIt;

      for(mapIt = this->mBcs.begin(); mapIt != this->mBcs.end(); ++mapIt)
      {
         tagMap.insert(std::make_pair(IoTools::IdToHuman::toTag(mapIt->first),mapIt->second));
      }

      return tagMap;
   }

}
