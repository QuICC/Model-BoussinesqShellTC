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
      : IntegratorBranch2D(physId, intgPhys, intgSpec, specId, fieldId, arithId), mIntgPart(intgPart)
   {
   }

   IntegratorBranch3D::~IntegratorBranch3D()
   {
   }

   IntgPartId IntegratorBranch3D::intgPartId() const
   {
      return this->mIntgPart;
   }

}
}
