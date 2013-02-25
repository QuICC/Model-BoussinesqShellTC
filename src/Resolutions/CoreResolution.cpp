/** \file CoreResolution.cpp
 *  \brief Source of the resolution object for a single CPU
 */

// System includes
//
#include <assert.h>

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

   SharedCTransformResolution CoreResolution::dim(const int i) const
   {
      // Check for correct sizes
      assert(i < this->mTransforms.size());

      return this->mTransforms.at(i);
   }

}
