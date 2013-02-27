/** \file SimulationResolution.cpp
 *  \brief Source of the simulation resolution object
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "Resolutions/SimulationResolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationResolution::SimulationResolution(const ArrayI& phys, const ArrayI& spec)
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

}
