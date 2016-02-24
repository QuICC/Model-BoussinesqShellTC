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
      : ProjectorBranch2D(specId, projSpec, projPhys, physId, fieldId, arithId), mProjPart(projPart)
   {
   }

   ProjectorBranch3D::ProjectorBranch3D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPartId projPart, ProjPhysId projPhys, std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id> physId, FieldType::Id fieldId, Arithmetics::Id arithId)
      : ProjectorBranch2D(specId, projSpec, projPhys, physId, fieldId, arithId), mProjPart(projPart)
   {
   }

   ProjectorBranch3D::~ProjectorBranch3D()
   {
   }

   ProjPartId ProjectorBranch3D::projPartId() const
   {
      return this->mProjPart;
   }

}
}
