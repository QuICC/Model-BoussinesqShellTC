/** 
 * @file SimulationResolution.cpp
 * @brief Source of the simulation resolution object
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Resolutions/SimulationResolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationResolution::SimulationResolution(const ArrayI& phys, const ArrayI& spec)
      : mBoxScale(Array::Ones(spec.size()))
   {
      // Assert both spaces have same dimensionalisty
      assert(phys.size() == spec.size());

      // Add physical space dimensions
      this->mDim.insert(std::make_pair(Dimensions::Space::PHYSICAL, phys));
      // Add spectral space dimensions
      this->mDim.insert(std::make_pair(Dimensions::Space::SPECTRAL, spec));
   }

   SimulationResolution::~SimulationResolution()
   {
   }

   int SimulationResolution::dim(const Dimensions::Simulation::Id simId, const Dimensions::Space::Id spaceId) const
   {
      // Safety assertion
      assert(this->mDim.find(spaceId)->second.size() > static_cast<int>(simId));

      return this->mDim.find(spaceId)->second(static_cast<int>(simId));
   }

   const ArrayI& SimulationResolution::dimensions(const Dimensions::Space::Id spaceId) const
   {
      return this->mDim.find(spaceId)->second;
   }

   MHDFloat SimulationResolution::boxScale(const Dimensions::Simulation::Id id) const
   {
      return this->mBoxScale(static_cast<Dimensions::Simulation::Id>(id));
   }

   void SimulationResolution::setBoxScale(const Array& boxScale)
   {
      this->mBoxScale = boxScale;
   }

}
