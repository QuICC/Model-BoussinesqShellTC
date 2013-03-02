/** \file FieldRequirement.cpp
 *  \brief Source of the variable requirements
 */

// System includes
//

// External includes
//

// Class include
//
#include "Variables/FieldRequirement.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   FieldRequirement::FieldRequirement(const bool isScalar, const bool needSpectral, const bool needPhysical, const bool needDiff)
      : mIsScalar(isScalar), mNeedSpectral(needSpectral), mNeedPhysical(needPhysical), mNeedDiff(needDiff)
   {
   }

   FieldRequirement::~FieldRequirement()
   {
   }

   bool FieldRequirement::isScalar() const
   {
      return this->mIsScalar;
   }

   bool FieldRequirement::needSpectral() const
   {
      return this->mNeedSpectral;
   }

   bool FieldRequirement::needPhysical() const
   {
      return this->mNeedPhysical;
   }

   bool FieldRequirement::needPhysicalDiff() const
   {
      return this->mNeedDiff;
   }

}
