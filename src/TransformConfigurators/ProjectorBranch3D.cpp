/** 
 * @file ProjectorBranch3D.cpp
 * @brief Source of the implementation of the backward tranform tree branch in 3D space
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
#include "TransformConfigurators/ProjectorBranch3D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorBranch3D::ProjectorBranch3D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPartId projPart, ProjPhysId projPhys, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mSpecId(specId), mProjSpec(projSpec), mProjPart(projPart), mProjPhys(projPhys), mPhysId(physId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   ProjectorBranch3D::~ProjectorBranch3D()
   {
   }

   FieldComponents::Spectral::Id ProjectorBranch3D::specId() const
   {
      return this->mSpecId;
   }

   ProjSpecId ProjectorBranch3D::projSpecId() const
   {
      return this->mProjSpec;
   }

   ProjPartId ProjectorBranch3D::projPartId() const
   {
      return this->mProjPart;
   }

   ProjPhysId ProjectorBranch3D::projPhysId() const
   {
      return this->mProjPhys;
   }

   FieldComponents::Physical::Id ProjectorBranch3D::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id ProjectorBranch3D::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id ProjectorBranch3D::arithId() const
   {
      return this->mArithId;
   }

}
}
