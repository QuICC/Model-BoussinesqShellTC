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

// Project includes
//
#include "TypeSelectors/TransformLeafSelector.hpp"
#include "TransformConfigurators/IntegratorBranch2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class describes a forward transform tree branch for 3D space
    */
   class IntegratorBranch3D: public IntegratorBranch2D
   {  
      public:
         /**
          * @brief Contructor for branch
          */
         IntegratorBranch3D(FieldComponents::Physical::Id physId, IntgPhysId intgPhys, IntgPartId intgPart, IntgSpecId intgSpec, FieldComponents::Spectral::Id specId, FieldType::Id fieldId, Arithmetics::Id arithId = Arithmetics::SET);

         /**
          * @brief Destructor
          */
         ~IntegratorBranch3D();

         /**
          * @brief Get 2D transform operator ID
          */
         IntgPartId intgPartId() const;

      private:
         /**
          * @brief Intgection operation in second dimension
          */
         IntgPartId  mIntgPart;
   };

}
}

#endif // INTEGRATORBRANCH3D_HPP
