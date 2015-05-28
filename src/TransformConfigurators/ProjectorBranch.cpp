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

   ProjectorBranch::ProjectorBranch(FieldComponents::Spectral::Id specId, ProjectorBranch::Proj1DId proj1D, ProjectorBranch::Proj2DId proj2D, ProjectorBranch::Proj3DId proj3D, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId)
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

   ProjectorBranch::Proj1DId ProjectorBranch::proj1DId() const
   {
      return this->mProj1D;
   }

   ProjectorBranch::Proj2DId ProjectorBranch::proj2DId() const
   {
      return this->mProj2D;
   }

   ProjectorBranch::Proj3DId ProjectorBranch::proj3DId() const
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
