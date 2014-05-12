/**
 * @file FieldIds.hpp
 * @brief Definition of some useful enums used to access fields by ID 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FIELDIDS_HPP
#define FIELDIDS_HPP

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
         DENSITY,
         /// Codensity field
         CODENSITY,
         /// Pressure field
         PRESSURE,
         /// Temperature field
         TEMPERATURE,
         /// Mean temperature field
         MEANTEMPERATURE,
         /// Streamfunction field
         STREAMFUNCTION,
         /// Axial velocity field
         VELOCITYZ,
         /// Axial vorticity field
         VORTICITYZ,
         /// Magnetic field
         MAGNETIC,
         /// Velocity field
         VELOCITY, 
         /// Vorticity field
         VORTICITY
      };
   };

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
         enum Id {
            /// First vector component
            ONE,
            /// Second vector component
            TWO,
            /// Third vector component
            THREE,
            /// Is a scalar
            SCALAR,
            /// Is not used
            NOTUSED
         };
      };

      /**
       * @brief Struct for the spectral field components
       */
      struct Spectral
      {
         /**
          * @brief Enum for Spectral field vector components
          */
         enum Id {
            /// First spectral vector component
            ONE,
            /// Second spectral vector component
            TWO,
            /// Third spectral vector component
            THREE,
            /// Is spectral scalar
            SCALAR,
            /// Is not used
            NOTUSED
         };
      };
   };

   /// Typedef for a full ID for a spectral field component
   typedef std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id>   SpectralFieldId;

   /// Typedef for a full ID for a physical field component
   typedef std::pair<PhysicalNames::Id,FieldComponents::Physical::Id>   PhysicalFieldId;
}

#endif // FIELDIDS_HPP
