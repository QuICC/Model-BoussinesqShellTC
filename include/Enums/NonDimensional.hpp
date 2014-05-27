/**
 * @file NonDimensional.hpp
 * @brief Definition of some useful enums for nondimensional paramters 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NONDIMENSIONAL_HPP
#define NONDIMENSIONAL_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

      /**
       * @brief Simple struct holding the mapping for the nondimensional parameters
       */
      struct NonDimensional {
         /**
          * @brief Enums of the different nondimensional factors
          */
         enum Id {
            /// Chi
            CHI,
            /// Ekman number
            EKMAN,
            /// Gamma
            GAMMA,
            /// Gap width R_o - R_i
            GAPWIDTH,
            /// Magnetic Ekman number
            MAGEKMAN,
            /// Magnetic Prandtl number
            MAGPRANDTL,
            /// Prandtl number
            PRANDTL,
            /// Rayleigh number
            RAYLEIGH,
            /// Roberts number
            ROBERTS,
            /// Rossby number
            ROSSBY,
            /// Radii ratio R_i/R_o
            RRATIO,
            /// Taylor number
            TAYLOR,
            /// Theta
            THETA,
         };
      };
}

#endif // NONDIMENSIONAL_HPP
