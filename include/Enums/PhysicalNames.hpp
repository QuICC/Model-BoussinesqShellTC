/** \file PhysicalNames.hpp
 *  \brief Definition of some useful enums mapping physical names to ID
 *
 *  \mhdBug Needs test
 */

#ifndef PHYSICALNAMES_HPP
#define PHYSICALNAMES_HPP

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
    * @brief Simple struct to hold the physical names IDs
    */
   struct PhysicalNames
   {
      /**
       * @name Enum for physical name to ID mapping
       */
      enum Id {CODENSITY, PRESSURE, TEMPERATURE, STREAMFUNCTION, VELOCITYZ, VORTICITYZ, PHI, MAGNETIC, VELOCITY, VORTICITY};
   };
}

#endif // PHYSICALNAMES_HPP
