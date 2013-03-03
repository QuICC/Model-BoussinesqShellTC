/** \file IdToHuman.hpp
 *  \brief Small routines to convert enum ID into human strings
 */

#ifndef IDTOHUMAN_HPP
#define IDTOHUMAN_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.hpp"

namespace GeoMHDiSCC {

namespace IoTools {

   /**
    * @brief Simple struct to hold the physical names IDs
    */
   struct IdToHuman
   {
      /**
       * @brief Convert ID to string
       */
      static std::string toString(const PhysicalNames::Id id);

      /**
       * @brief Convert ID to string
       */
      static std::string toString(const FieldComponents::Physical::Id id);

      /**
       * @brief Convert ID to string tag
       */
      static std::string toTag(const FieldComponents::Physical::Id id);

      /**
       * @brief Convert ID to string
       */
      static std::string toString(const FieldComponents::Spectral::Id id);

      /**
       * @brief Convert ID to string tag
       */
      static std::string toTag(const FieldComponents::Spectral::Id id);

      private:
         /**
          * @brief Constructor
          */
         IdToHuman();

         /**
          * @brief Constructor
          */
         ~IdToHuman();
   };
}
}

#endif // IDTOHUMAN_HPP
