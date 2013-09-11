/** 
 * @file IdToHuman.hpp
 * @brief Small routines to convert enum ID into human strings
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"

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
       * @brief Convert ID to tag
       */
      static std::string toTag(const PhysicalNames::Id id);

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

      /**
       * @brief Convert ID to tag
       */
      static std::string toTag(const NonDimensional::Id id);

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
