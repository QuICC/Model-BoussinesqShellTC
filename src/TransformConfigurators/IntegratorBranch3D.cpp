/** 
 * @file IntegratorBranch3D.cpp
 * @brief Source of the implementation of the forward transform tree branch for 3D
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
#include "TransformConfigurators/IntegratorBranch3D.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorBranch3D::IntegratorBranch3D(FieldComponents::Physical::Id physId, IntgPhysId intgPhys, IntgPartId intgPart, IntgSpecId intgSpec, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mPhysId(physId), mIntgSpec(intgSpec), mIntgPart(intgPart), mIntgPhys(intgPhys), mSpecId(specId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   IntegratorBranch3D::~IntegratorBranch3D()
   {
   }

   FieldComponents::Spectral::Id IntegratorBranch3D::specId() const
   {
      return this->mSpecId;
   }

   IntgSpecId IntegratorBranch3D::intgSpecId() const
   {
      return this->mIntgSpec;
   }

   IntgPartId IntegratorBranch3D::intgPartId() const
   {
      return this->mIntgPart;
   }

   IntgPhysId IntegratorBranch3D::intgPhysId() const
   {
      return this->mIntgPhys;
   }

   FieldComponents::Physical::Id IntegratorBranch3D::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id IntegratorBranch3D::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id IntegratorBranch3D::arithId() const
   {
      return this->mArithId;
   }

}
}
