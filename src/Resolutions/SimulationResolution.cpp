/** 
 * @file SimulationResolution.cpp
 * @brief Source of the simulation resolution object
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

   SimulationResolution::SimulationResolution(const ArrayI& phys, const ArrayI& spec, const ArrayI& trans)
      : mBoxScale(Array::Ones(spec.size()))
   {
      // Assert both spaces have same dimensionalisty
      assert(phys.size() == spec.size());

      // Add physical space dimensions
      this->mDim.insert(std::make_pair(Dimensions::Space::PHYSICAL, phys));
      // Add spectral space dimensions
      this->mDim.insert(std::make_pair(Dimensions::Space::SPECTRAL, spec));
      // Add transform space dimensions
      this->mDim.insert(std::make_pair(Dimensions::Space::TRANSFORM, trans));
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
