/** 
 * @file CoreResolution.cpp
 * @brief Source of the resolution object for a single CPU
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
      assert(static_cast<size_t>(id) < this->mTransforms.size());

      return this->mTransforms.at(static_cast<int>(id));
   }

}
