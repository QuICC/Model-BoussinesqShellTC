/** \file CoreResolution.cpp
 *  \brief Source of the resolution object for a single CPU
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "Resolutions/CoreResolution.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   CoreResolution::CoreResolution(const std::vector<SharedTransformResolution>& transformRes)
      : mTransforms(transformRes)
   {
   }

   CoreResolution::~CoreResolution()
   {
   }

   int CoreResolution::nDim() const
   {
      return this->mTransforms.size();
   }

   SharedCTransformResolution CoreResolution::dim(const Dimensions::Transform::Id id) const
   {
      // Check for correct sizes
      assert(static_cast<int>(id) >= 0);
      assert(static_cast<unsigned int>(id) < this->mTransforms.size());

      return this->mTransforms.at(static_cast<int>(id));
   }

}
