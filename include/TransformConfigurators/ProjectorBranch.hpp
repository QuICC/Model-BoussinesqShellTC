/** 
 * @file ProjectorBranch.hpp
 * @brief This class defines a backward transform branch
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PROJECTORBRANCH_HPP
#define PROJECTORBRANCH_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TreeSelector.hpp"
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class ProjectorBranch
   {  
      public:
         /**
          * @brief Contructor for branch
          */
         ProjectorBranch(FieldComponents::Spectral::Id specId, ProjSpecId proj1D, ProjPartId proj2D, ProjPhysId proj3D, FieldComponents::Physical::Id physId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Destructor
          */
         ~ProjectorBranch();

         /**
          * @brief Get spectral component ID
          */
         FieldComponents::Spectral::Id specId() const;

         /**
          * @brief Get 1D transform operator ID
          */
         ProjSpecId proj1DId() const;

         /**
          * @brief Get 2D transform operator ID
          */
         ProjPartId proj2DId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         ProjPhysId proj3DId() const;

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
         ProjSpecId  mProj1D;

         /**
          * @brief Projection operation in fourth dimension
          */
         ProjPartId  mProj2D;

         /**
          * @brief Projection operation in third dimension
          */
         ProjPhysId  mProj3D;

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

#endif // PROJECTORBRANCH_HPP
