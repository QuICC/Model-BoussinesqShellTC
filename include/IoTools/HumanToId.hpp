/** 
 * @file HumanToId.hpp
 * @brief Small routines to convert human string to enum ID
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef HUMANTOID_HPP
#define HUMANTOID_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/ModelOperator.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   /**
    *  @brief Small routines to convert human string to enum ID
    */
   struct HumanToId
   {
      /**
       * @brief Convert tag to nondimensional parameter ID
       */
      static NonDimensional::Id toNd(const std::string& id);

      /**
       * @brief Convert string to model operator ID
       */
      static ModelOperator::Id toModOp(const std::string& id);

      /**
       * @brief Convert string to physical name ID
       */
      static PhysicalNames::Id toPhys(const std::string& id);

      /**
       * @brief Convert string to field component ID
       */
      static FieldComponents::Spectral::Id toComp(const std::string& id);

      private:
         /**
          * @brief Constructor
          */
         HumanToId();

         /**
          * @brief Constructor
          */
         ~HumanToId();
   };
}
}

#endif // HUMANTOID_HPP
