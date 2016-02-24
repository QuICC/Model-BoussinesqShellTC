/** 
 * @file ProjectorBranch2D.cpp
 * @brief Source of the implementation of the backward tranform tree branch in 2D space
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
#include "TransformConfigurators/ProjectorBranch2D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   ProjectorBranch2D::ProjectorBranch2D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPhysId projPhys, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mSpecId(specId), mProjSpec(projSpec), mProjPhys(projPhys), mPhysId(physId), mFieldId(fieldId), mArithId(arithId)
   {
   }

   ProjectorBranch2D::ProjectorBranch2D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPhysId projPhys, std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id> physId, FieldType::Id fieldId, Arithmetics::Id arithId)
      :mSpecId(specId), mProjSpec(projSpec), mProjPhys(projPhys), mPhysId(physId.first), mFieldId(fieldId), mArithId(arithId)
   {
   }

   ProjectorBranch2D::~ProjectorBranch2D()
   {
   }

   FieldComponents::Spectral::Id ProjectorBranch2D::specId() const
   {
      return this->mSpecId;
   }

   ProjSpecId ProjectorBranch2D::projSpecId() const
   {
      return this->mProjSpec;
   }

   ProjPhysId ProjectorBranch2D::projPhysId() const
   {
      return this->mProjPhys;
   }

   FieldComponents::Physical::Id ProjectorBranch2D::physId() const
   {
      return this->mPhysId;
   }

   FieldType::Id ProjectorBranch2D::fieldId() const
   {
      return this->mFieldId;
   }

   Arithmetics::Id ProjectorBranch2D::arithId() const
   {
      return this->mArithId;
   }

}
}
