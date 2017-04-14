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

namespace QuICC {

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
         /// large scale x B field 
         BX,
         /// large scale y B field 
         BY,
         /// Density field
         DENSITY,
         /// X derivative of mean temperature field
         DX_MEANTEMPERATURE,
         /// Z derivative of mean temperature field
         DZ_MEANTEMPERATURE,
         /// Entropy field
         ENTROPY,
         /// Electromagnetic force x (in By equation)
         EMFX,
         /// Electromagnetic force y (in Bx equation)
         EMFY,
         /// Fluctuating bx
         FBX,
         /// Fluctuating by
         FBY,
         /// Fluctuating bz
         FBZ,
         /// Magnetic field
         FJZ,
         /// Z fluctuating density current
         MAGNETIC,
         /// X Magnetic field
         MAGNETICX,
         /// Y Magnetic field
         MAGNETICY,
         /// Z Magnetic field
         MAGNETICZ,
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
         /// Zonal velocity field (m = 0)
         ZONAL_VELOCITY, 
         /// Non-zonal velocity field (m != 0)
         NONZONAL_VELOCITY, 
         /// Vorticity field
         VORTICITY,
         /// X vorticity field
         VORTICITYX,
         /// Y vorticity field
         VORTICITYY,
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

         /// Fluctuating temperature field
         FLUCT_TEMPERATURE,
         /// Fluctuating X magnetic field
         FLUCT_MAGNETIC,
         /// Fluctuating X magnetic field
         FLUCT_MAGNETICX,
         /// Fluctuating Y magnetic field
         FLUCT_MAGNETICY,
         /// Fluctuating Z magnetic field
         FLUCT_MAGNETICZ,
         /// Fluctuating X velocity field
         FLUCT_VELOCITY,
         /// Fluctuating X velocity field
         FLUCT_VELOCITYX,
         /// Fluctuating Y velocity field
         FLUCT_VELOCITYY,
         /// Fluctuating Z velocity field
         FLUCT_VELOCITYZ,

         /// Mean temperature field
         MEAN_TEMPERATURE,
         /// Mean X magnetic field
         MEAN_MAGNETIC,
         /// Mean X magnetic field
         MEAN_MAGNETICX,
         /// Mean Y magnetic field
         MEAN_MAGNETICY,
         /// Mean Z magnetic field
         MEAN_MAGNETICZ,
         /// Mean X velocity field
         MEAN_VELOCITY,
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

	      /// Imposed magnetic field
	      IMPOSED_MAGNETIC,
      };
   };

   /**
    * @brief Simple struct to hold the field types
    */
   struct FieldType
   {
      /**
       * @brief Enum for the field types
       */
      enum Id {
         /// Scalar field
         SCALAR,
         /// Vector field
         VECTOR,
         /// Gradient field
         GRADIENT,
         /// Curl field
         CURL,
         /// 2nd order gradient field
         GRADIENT2,
         /// divergence field
         DIVERGENCE,
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
            /// X component of cartesian field
            X,
            /// Radial component of cylindrical or spherical field
            R,

            /// Y component of cartesian field
            Y,
            /// Theta component of cylindrical or spherical field
            THETA,

            /// Z component of cartesian field
            Z,
            /// Phi component of spherical field
            PHI,

            /// Is a scalar
            SCALAR,

            /// Is not used
            NOTUSED,

            // Define generic enums for cartesian geometry
            #if defined QUICC_SPATIALSCHEME_TTT || defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF 
            /// First vector component
            ONE = X,
            /// Second vector component
            TWO = Y,
            /// Third vector component
            THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_TFF
            /// First vector component
            ONE = Z,
            /// Second vector component
            TWO = X,
            /// Third vector component
            THREE = Y,

            #elif defined QUICC_SPATIALSCHEME_WFT
            /// First vector component
            ONE = R,
            /// Second vector component
            TWO = THETA,
            /// Third vector component
            THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_CFT || defined QUICC_SPATIALSCHEME_AFT
            /// First vector component
            ONE = R,
            /// Second vector component
            TWO = THETA,
            /// Third vector component
            THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFL || defined QUICC_SPATIALSCHEME_WLFM
            /// First vector component
            ONE = R,
            /// Second vector component
            TWO = THETA,
            /// Third vector component
            THREE = PHI,

            #elif defined QUICC_SPATIALSCHEME_TT
            /// First vector component
            ONE = X,
            /// Second vector component
            TWO = Z,
            /// Third vector component
            THREE = NOTUSED,

            #elif defined QUICC_SPATIALSCHEME_TF
            /// First vector component
            ONE = Z,
            /// Second vector component
            TWO = X,
            /// Third vector component
            THREE = NOTUSED,

            #elif defined QUICC_SPATIALSCHEME_CF || defined QUICC_SPATIALSCHEME_AF
            /// First vector component
            ONE = R,
            /// Second vector component
            TWO = THETA,
            /// Third vector component
            THREE = NOTUSED,

            #endif // defined QUICC_SPATIALSCHEME_TTT || defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF
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
            /// X component of cartesian field
            X,
            /// Radial component of cylindrical or spherical field
            R,
            /// Toroidal component
            TOR,
            /// Q component of Spheroidal/Toroidal
            Q,

            /// Y component of cartesian field
            Y,
            /// Theta component of cylindrical or spherical field
            THETA,
            /// Poloidal component of Toroidal/Poloidal
            POL,
            /// S component of Spheroidal/Toroidal
            S,

            /// Z component of cartesian field
            Z,
            /// Phi component of spherical field
            PHI,
            /// T component of spheroidal/toroidal field
            T,

            /// Is spectral scalar
            SCALAR,

            /// Is not used
            NOTUSED,

            // Define generic enums for cartesian geometry
            #if defined QUICC_SPATIALSCHEME_TFF_TORPOL || defined QUICC_SPATIALSCHEME_SLFL_TORPOL || defined QUICC_SPATIALSCHEME_SLFM_TORPOL || defined QUICC_SPATIALSCHEME_BLFL_TORPOL || defined QUICC_SPATIALSCHEME_BLFM_TORPOL || defined QUICC_SPATIALSCHEME_WLFL_TORPOL || defined QUICC_SPATIALSCHEME_WLFM_TORPOL || defined QUICC_SPATIALSCHEME_WFT_TORPOL
               /// First vector component
               ONE = TOR,
               /// Second vector component
               TWO = POL,
               /// Third vector component
               THREE = NOTUSED,

            #elif defined QUICC_SPATIALSCHEME_SLFL_QST || defined QUICC_SPATIALSCHEME_SLFM_QST || defined QUICC_SPATIALSCHEME_BLFL_QST || defined QUICC_SPATIALSCHEME_BLFM_QST || defined QUICC_SPATIALSCHEME_WLFL_QST || defined QUICC_SPATIALSCHEME_WLFM_QST
               /// First vector component
               ONE = Q,
               /// Second vector component
               TWO = S,
               /// Third vector component
               THREE = T,

            #elif defined QUICC_SPATIALSCHEME_TTT || defined QUICC_SPATIALSCHEME_TFT || defined QUICC_SPATIALSCHEME_FFF 
               /// First vector component
               ONE = X,
               /// Second vector component
               TWO = Y,
               /// Third vector component
               THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_TFF
               /// First vector component
               ONE = Z,
               /// Second vector component
               TWO = X,
               /// Third vector component
               THREE = Y,

            #elif defined QUICC_SPATIALSCHEME_WFT
               /// First vector component
               ONE = R,
               /// Second vector component
               TWO = THETA,
               /// Third vector component
               THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_CFT || defined QUICC_SPATIALSCHEME_AFT
               /// First vector component
               ONE = R,
               /// Second vector component
               TWO = THETA,
               /// Third vector component
               THREE = Z,

            #elif defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFL || defined QUICC_SPATIALSCHEME_WLFM
               /// First vector component
               ONE = R,
               /// Second vector component
               TWO = THETA,
               /// Third vector component
               THREE = PHI,

            #elif defined QUICC_SPATIALSCHEME_TT
               /// First vector component
               ONE = X,
               /// Second vector component
               TWO = Z,
               /// Third vector component
               THREE = NOTUSED,

            #elif defined QUICC_SPATIALSCHEME_TF
               /// First vector component
               ONE = Z,
               /// Second vector component
               TWO = X,
               /// Third vector component
               THREE = NOTUSED,

            #elif defined QUICC_SPATIALSCHEME_AF || defined QUICC_SPATIALSCHEME_CF 
               /// First vector component
               ONE = R,
               /// Second vector component
               TWO = THETA,
               /// Third vector component
               THREE = NOTUSED,
            #endif //defined QUICC_SPATIALSCHEME_TFF_TORPOL || defined QUICC_SPATIALSCHEME_SLFL_TORPOL || defined QUICC_SPATIALSCHEME_SLFM_TORPOL || defined QUICC_SPATIALSCHEME_BLFL_TORPOL || defined QUICC_SPATIALSCHEME_BLFM_TORPOL || defined QUICC_SPATIALSCHEME_WLFL_TORPOL || defined QUICC_SPATIALSCHEME_WLFM_TORPOL || defined QUICC_SPATIALSCHEME_WFT_TORPOL
         };
      };
   };

   /// Typedef for a full ID for a spectral field component
   typedef std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id>   SpectralFieldId;

   /// Typedef for a full ID for a physical field component
   typedef std::pair<PhysicalNames::Id,FieldComponents::Physical::Id>   PhysicalFieldId;
}

#endif // FIELDIDS_HPP
