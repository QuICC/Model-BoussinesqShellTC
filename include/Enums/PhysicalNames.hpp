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
      enum Id {
         /// Codensity field
         CODENSITY,
         /// Pressure field
         PRESSURE,
         /// Temperature field
         TEMPERATURE,
         /// Streamfunction field
         STREAMFUNCTION,
         /// Axial velocity field
         VELOCITYZ,
         /// Axial vorticity field
         VORTICITYZ,
         /// Phi field
         PHI,
         /// Magnetic field
         MAGNETIC,
         /// Velocity field
         VELOCITY, 
         /// Vorticity field
         VORTICITY
      };
   };
}

#endif // PHYSICALNAMES_HPP
