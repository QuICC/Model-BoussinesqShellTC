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
#include "TypeSelectors/TreeSelector.hpp"
#include "Enums/Arithmetics.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a forward transform tree branch
    */
   class IntegratorBranch
   {  
      public:
         /**
          * @brief Contructor for branch
          */
         IntegratorBranch(FieldComponents::Physical::Id physId, IntgPhysId intg3D, IntgPartId intg2D, IntgSpecId intg1D, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

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
         IntgSpecId intg1DId() const;

         /**
          * @brief Get 2D transform operator ID
          */
         IntgPartId intg2DId() const;

         /**
          * @brief Get 3D transform operator ID
          */
         IntgPhysId intg3DId() const;

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
         IntgSpecId  mIntg1D;

         /**
          * @brief Intgection operation in fourth dimension
          */
         IntgPartId  mIntg2D;

         /**
          * @brief Intgection operation in third dimension
          */
         IntgPhysId  mIntg3D;
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
