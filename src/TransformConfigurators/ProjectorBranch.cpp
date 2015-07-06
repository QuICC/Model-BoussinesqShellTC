/** 
 * @file ProjectorBranch.cpp
 * @brief Source of the implementation of the backward tranform tree branch
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
#include "TransformConfigurators/ProjectorBranch.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorBranch::ProjectorBranch(FieldComponents::Spectral::Id specId, ProjSpecId proj1D, ProjPartId proj2D, ProjPhysId proj3D, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mSpecId(specId), mProj1D(proj1D), mProj2D(proj2D), mProj3D(proj3D), mPhysId(physId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   ProjectorBranch::~ProjectorBranch()
   {
   }

   FieldComponents::Spectral::Id ProjectorBranch::specId() const
   {
      return this->mSpecId;
   }

   ProjSpecId ProjectorBranch::proj1DId() const
   {
      return this->mProj1D;
   }

   ProjPartId ProjectorBranch::proj2DId() const
   {
      return this->mProj2D;
   }

   ProjPhysId ProjectorBranch::proj3DId() const
   {
      return this->mProj3D;
   }

   FieldComponents::Physical::Id ProjectorBranch::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id ProjectorBranch::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id ProjectorBranch::arithId() const
   {
      return this->mArithId;
   }

}
}
