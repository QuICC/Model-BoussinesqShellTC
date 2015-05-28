/** 
 * @file IntegratorBranch.cpp
 * @brief Source of the implementation of the forward transform tree branch
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/IntegratorBranch.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorBranch::IntegratorBranch(FieldComponents::Physical::Id physId, IntegratorBranch::Intg3DId intg3D, IntegratorBranch::Intg2DId intg2D, IntegratorBranch::Intg1DId intg1D, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mPhysId(physId), mIntg1D(intg1D), mIntg2D(intg2D), mIntg3D(intg3D), mSpecId(specId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   IntegratorBranch::~IntegratorBranch()
   {
   }

   FieldComponents::Spectral::Id IntegratorBranch::specId() const
   {
      return this->mSpecId;
   }

   IntegratorBranch::Intg1DId IntegratorBranch::intg1DId() const
   {
      return this->mIntg1D;
   }

   IntegratorBranch::Intg2DId IntegratorBranch::intg2DId() const
   {
      return this->mIntg2D;
   }

   IntegratorBranch::Intg3DId IntegratorBranch::intg3DId() const
   {
      return this->mIntg3D;
   }

   FieldComponents::Physical::Id IntegratorBranch::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id IntegratorBranch::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id IntegratorBranch::arithId() const
   {
      return this->mArithId;
   }

}
}
