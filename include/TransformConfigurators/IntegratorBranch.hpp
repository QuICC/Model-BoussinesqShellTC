/** 
 * @file IntegratorBranch.hpp
 * @brief This class defines a forward transform branch
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORBRANCH_HPP
#define INTEGRATORBRANCH_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Intgect includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a forward transform tree branch
    */
   class IntegratorBranch
   {  
      public:
         /// Typedefs to simplify definition of intgection operators
         typedef TransformCoordinatorType::Transform1DType::IntegratorType Intg1DType;
         typedef TransformCoordinatorType::Transform2DType::IntegratorType Intg2DType;
         typedef TransformCoordinatorType::Transform3DType::IntegratorType Intg3DType;
         typedef TransformCoordinatorType::Transform1DType::IntegratorType::Id Intg1DId;
         typedef TransformCoordinatorType::Transform2DType::IntegratorType::Id Intg2DId;
         typedef TransformCoordinatorType::Transform3DType::IntegratorType::Id Intg3DId;

         /**
          * @brief Contructor for branch
          */
         IntegratorBranch(FieldComponents::Physical::Id physId, Intg3DId intg3D, Intg2DId intg2D, Intg1DId intg1D, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Get physical component ID
          */
         FieldComponents::Physical::Id physId() const;

         /**
          * @brief Destructor
          */
         ~IntegratorBranch();

         /**
          * @brief Get 1D transform operator ID
          */
         Intg1DId intg1DId() const;

         /**
          * @brief Get 2D transform operator ID
          */
         Intg2DId intg2DId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         Intg3DId intg3DId() const;

         /**
          * @brief Get spectral component ID
          */
         FieldComponents::Spectral::Id specId() const;

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
          * @brief Physical component required for transform branch
          */
         FieldComponents::Physical::Id mPhysId;

         /**
          * @brief Intgection operation in first dimension
          */
         Intg1DId  mIntg1D;

         /**
          * @brief Intgection operation in fourth dimension
          */
         Intg2DId  mIntg2D;

         /**
          * @brief Intgection operation in third dimension
          */
         Intg3DId  mIntg3D;
         /**
          * @brief Spectral component required for transform branch
          */
         FieldComponents::Spectral::Id mSpecId;

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

#endif // INTEGRATORBRANCH_HPP
