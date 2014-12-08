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
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class ProjectorBranch
   {  
      public:
         /// Typedefs to simplify definition of projection operators
         typedef TransformCoordinatorType::Transform1DType::ProjectorType Proj1DType;
         typedef TransformCoordinatorType::Transform2DType::ProjectorType Proj2DType;
         typedef TransformCoordinatorType::Transform3DType::ProjectorType Proj3DType;
         typedef TransformCoordinatorType::Transform1DType::ProjectorType::Id Proj1DId;
         typedef TransformCoordinatorType::Transform2DType::ProjectorType::Id Proj2DId;
         typedef TransformCoordinatorType::Transform3DType::ProjectorType::Id Proj3DId;

         /**
          * @brief Contructor for branch
          */
         ProjectorBranch(FieldComponents::Spectral::Id specId, Proj1DId proj1D, Proj2DId proj2D, Proj3DId proj3D, FieldComponents::Physical::Id physId, FieldType::Id fieldId);

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
         Proj1DId proj1DId() const;

         /**
          * @brief Get 2D transform operator ID
          */
         Proj2DId proj2DId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         Proj3DId proj3DId() const;

         /**
          * @brief Get physical component ID
          */
         FieldComponents::Physical::Id physId() const;

         /**
          * @brief Get field type ID
          */
         FieldType::Id fieldId() const;

      private:
         /**
          * @brief Spectral component required for transform branch
          */
         FieldComponents::Spectral::Id mSpecId;

         /**
          * @brief Projection operation in first dimension
          */
         Proj1DId  mProj1D;

         /**
          * @brief Projection operation in fourth dimension
          */
         Proj2DId  mProj2D;

         /**
          * @brief Projection operation in third dimension
          */
         Proj3DId  mProj3D;

         /**
          * @brief Physical component required for transform branch
          */
         FieldComponents::Physical::Id mPhysId;

         /**
          * @brief Field type required for transform branch
          */
         FieldType::Id mFieldId;
   };

}
}

#endif // PROJECTORBRANCH_HPP
