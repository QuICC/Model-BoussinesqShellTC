/** 
 * @file ProjectorBranch2D.hpp
 * @brief This class defines a backward transform branch for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORBRANCH2D_HPP
#define PROJECTORBRANCH2D_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class ProjectorBranch2D
   {  
      public:
         /**
          * @brief Contructor for branch for single physical ID
          */
         ProjectorBranch2D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPhysId projPhys, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Contructor for branch for pair physical ID
          */
         ProjectorBranch2D(FieldComponents::Spectral::Id specId, ProjSpecId projSpec, ProjPhysId projPhys, std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id> physId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Destructor
          */
         ~ProjectorBranch2D();

         /**
          * @brief Get spectral component ID
          */
         FieldComponents::Spectral::Id specId() const;

         /**
          * @brief Get 1D transform operator ID
          */
         ProjSpecId projSpecId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         ProjPhysId projPhysId() const;

         /**
          * @brief Get physical component ID
          */
         FieldComponents::Physical::Id physId() const;

         /**
          * @brief Get field type ID
          */
         FieldType::Id fieldId() const;

         /**
          * @brief Get arithmetic ID
          */
         Arithmetics::Id arithId() const;

      private:
         /**
          * @brief Spectral component required for transform branch
          */
         FieldComponents::Spectral::Id mSpecId;

         /**
          * @brief Projection operation in first dimension
          */
         ProjSpecId  mProjSpec;

         /**
          * @brief Projection operation in third dimension
          */
         ProjPhysId  mProjPhys;

         /**
          * @brief Physical component required for transform branch
          */
         FieldComponents::Physical::Id mPhysId;

         /**
          * @brief Field type required for transform branch
          */
         FieldType::Id mFieldId;

         /**
          * @brief Arithmetic operation to store result after transform
          */
         Arithmetics::Id mArithId;
   };

}
}

#endif // PROJECTORBRANCH2D_HPP
