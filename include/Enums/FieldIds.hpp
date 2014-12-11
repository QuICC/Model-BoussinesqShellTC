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
         CODENSITY = 0,
         /// Density field
         DENSITY,
         /// X derivative of mean temperature field
         DX_MEANTEMPERATURE,
         /// Z derivative of mean temperature field
         DZ_MEANTEMPERATURE,
         /// Entropy field
         ENTROPY,
         /// Magnetic field
         MAGNETIC,
         /// Mean temperature field
         MEANTEMPERATURE,
         /// Pressure field
         PRESSURE,
         /// Temperature field
         TEMPERATURE,
         /// Streamfunction field
         STREAMFUNCTION,
         /// Velocity field
         VELOCITY, 
         /// X velocity field
         VELOCITYX,
         /// Y velocity field
         VELOCITYY,
         /// Z velocity field
         VELOCITYZ,
         /// Vorticity field
         VORTICITY,
         /// Axial vorticity field
         VORTICITYZ,

         /// Phi field (for example phi = D_z w)
         PHI,
   
         /// Non orthogonal streamfunction field
         NO_STREAMFUNCTION,
         /// Non orthogonal vertical velocity field
         NO_VELOCITYZ,
         /// Non orthogonal vertical vorticity field
         NO_VORTICITYZ,
   
         /// Temperature field in tilted box
         TILTED_TEMPERATURE,
         /// Vertical streamfunction field in tilted box
         TILTED_STREAMFUNCTION,
         /// vertical velocity field in tilted box
         TILTED_VELOCITYZ,
         /// Non orthogonal streamfunction field in tilted box
         TILTED_NO_STREAMFUNCTION,
         /// Non orthogonal vertical velocity field in tilted box
         TILTED_NO_VELOCITYZ,
         /// Non orthogonal vertical vorticity field in tilted box
         TILTED_NO_VORTICITYZ,

         /// Mean X velocity field
         MEAN_VELOCITYX,
         /// Mean Y velocity field
         MEAN_VELOCITYY,
         /// Mean Z velocity field
         MEAN_VELOCITYZ,

         /// kinetic energy: u \cdot u
         KINETIC_ENERGY,
         /// Zonal kinetic energy: u \cdot u
         ZONAL_KINETIC_ENERGY,
         /// Non zonal kinetic energy: u \cdot u
         NONZONAL_KINETIC_ENERGY,
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
