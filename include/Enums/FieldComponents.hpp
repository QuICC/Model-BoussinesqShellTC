/** \file FieldComponents.hpp
 *  \brief Definition of some useful enums used to access scalar fields and vector field components
 *
 *  \mhdBug Needs test
 */

#ifndef FIELDCOMPONENTS_HPP
#define FIELDCOMPONENTS_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

   /**
    * @brief Simple struct to hold the field components
    */
   struct FieldComponents
   {
      /**
       * @brief Struct for the physical field components
       */
      struct Physical
      {
         /**
          * @brief Enum for physical field vector components
          */
         enum Id {ONE, TWO, THREE, SCALAR, NOTUSED};
      };

      /**
       * @brief Struct for the spectral field components
       */
      struct Spectral
      {
         /**
          * @brief Enum for Spectral field vector components
          */
         enum Id {ONE, TWO, THREE, SCALAR, NOTUSED};
      };
   };
}

#endif // FIELDCOMPONENTS_HPP
