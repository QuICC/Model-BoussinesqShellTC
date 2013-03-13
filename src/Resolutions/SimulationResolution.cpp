/** \file SimulationResolution.cpp
 *  \brief Source of the simulation resolution object
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

   ArrayI SimulationResolution::orderedDims(const Dimensions::Space::Id spaceId) const
   {
      int nDims = this->mDim.find(spaceId)->second.size();
      ArrayI oDims(nDims);

      if(spaceId == Dimensions::Space::SPECTRAL)
      {
         // Get the dimension in 1D, 3D, 2D order
         oDims(0) = this->mDim.find(spaceId)->second(0);
         for(int i = 1; i < nDims; ++i)
         {
            oDims(i) = this->mDim.find(spaceId)->second(nDims-i);
         }
      } else //if(spaceId == Dimensions::Space::PHYSICAL)
      {
         // Get the dimensions in reverse order 3D, 2D, 1D
         for(int i = 0; i < nDims; ++i)
         {
            oDims(i) = this->mDim.find(spaceId)->second(nDims-1-i);
         }
      }

      return oDims;
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
