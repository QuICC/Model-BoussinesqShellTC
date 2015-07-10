/** 
 * @file IntegratorBranch2D.hpp
 * @brief This class defines a forward transform branch for 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INTEGRATORBRANCH2D_HPP
#define INTEGRATORBRANCH2D_HPP

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
    * @brief This class describes a forward transform tree branch for 3D space
    */
   class IntegratorBranch2D
   {  
      public:
         /**
          * @brief Contructor for branch
          */
         IntegratorBranch2D(FieldComponents::Physical::Id physId, IntgPhysId intgPhys, IntgSpecId intgSpec, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Get physical component ID
          */
         FieldComponents::Physical::Id physId() const;

         /**
          * @brief Destructor
          */
         ~IntegratorBranch2D();

         /**
          * @brief Get spectral transform operator ID
          */
         IntgSpecId intgSpecId() const;

         /**
          * @brief Get physical transform operator ID
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

#endif // INTEGRATORBRANCH2D_HPP
