/** 
 * @file IntegratorBranch3D.hpp
 * @brief This class defines a forward transform branch for 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORBRANCH3D_HPP
#define INTEGRATORBRANCH3D_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Intgect includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a forward transform tree branch for 3D space
    */
   class IntegratorBranch3D
   {  
      public:
         /**
          * @brief Contructor for branch
          */
         IntegratorBranch3D(FieldComponents::Physical::Id physId, IntgPhysId intgPhys, IntgPartId intgPart, IntgSpecId intgSpec, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Get physical component ID
          */
         FieldComponents::Physical::Id physId() const;

         /**
          * @brief Destructor
          */
         ~IntegratorBranch3D();

         /**
          * @brief Get 1D transform operator ID
          */
         IntgSpecId intgSpecId() const;

         /**
          * @brief Get 2D transform operator ID
          */
         IntgPartId intgPartId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         IntgPhysId intgPhysId() const;

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
         IntgSpecId  mIntgSpec;

         /**
          * @brief Intgection operation in second dimension
          */
         IntgPartId  mIntgPart;

         /**
          * @brief Intgection operation in third dimension
          */
         IntgPhysId  mIntgPhys;
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

#endif // INTEGRATORBRANCH3D_HPP
