/** 
 * @file IntegratorBranch2D.cpp
 * @brief Source of the implementation of the forward transform tree branch for 2D
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
#include "TransformConfigurators/IntegratorBranch2D.hpp"

// Intgect includes
//

namespace GeoMHDiSCC {

namespace Transform {

   IntegratorBranch2D::IntegratorBranch2D(FieldComponents::Physical::Id physId, IntgPhysId intgPhys, IntgSpecId intgSpec, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mPhysId(physId), mIntgSpec(intgSpec), mIntgPhys(intgPhys), mSpecId(specId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   IntegratorBranch2D::~IntegratorBranch2D()
   {
   }

   FieldComponents::Spectral::Id IntegratorBranch2D::specId() const
   {
      return this->mSpecId;
   }

   IntgSpecId IntegratorBranch2D::intgSpecId() const
   {
      return this->mIntgSpec;
   }

   IntgPhysId IntegratorBranch2D::intgPhysId() const
   {
      return this->mIntgPhys;
   }

   FieldComponents::Physical::Id IntegratorBranch2D::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id IntegratorBranch2D::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id IntegratorBranch2D::arithId() const
   {
      return this->mArithId;
   }

}
}
