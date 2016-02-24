/** 
 * @file ProjectorBranch3D.hpp
 * @brief This class defines a backward transform branch
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORBRANCH3D_HPP
#define PROJECTORBRANCH3D_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "TransformConfigurators/ProjectorBranch2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class ProjectorBranch3D: public ProjectorBranch2D
   {  
      public:
         /**
          * @brief Contructor for branch for single physical ID
          */
         ProjectorBranch3D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPartId projPart, ProjPhysId projPhys, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Contructor for branch for pair physical ID
          */
         ProjectorBranch3D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPartId projPart, ProjPhysId projPhys, std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id> physId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Destructor
          */
         ~ProjectorBranch3D();

         /**
          * @brief Get 2D transform operator ID
          */
         ProjPartId projPartId() const;

      private:
         /**
          * @brief Projection operation in second dimension
          */
         ProjPartId  mProjPart;
   };

}
}

#endif // PROJECTORBRANCH3D_HPP
