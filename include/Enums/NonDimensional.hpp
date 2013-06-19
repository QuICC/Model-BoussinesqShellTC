/** \file NonDimensional.hpp
 *  \brief Definition of some useful enums for nondimensional paramters
 */

#ifndef NONDIMESIONAL_HPP
#define NONDIMESIONAL_HPP

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
            /// Ekman number
            EKMAN,
            /// Roberts number
            ROBERTS,
            /// Rayleigh number
            RAYLEIGH,
            /// Rossby number
            ROSSBY,
            /// Magnetic Ekman number
            MAGEKMAN,
            /// Prandtl number
            PRANDTL,
            /// Magnetic Prandtl number
            MAGPRANDTL,
            /// Chi angle
            CHI,
            /// Topographic ratio
            GAMMA,
            /// Gap width R_o - R_i
            GAPWIDTH,
            /// Radii ratio R_i/R_o
            RRATIO,
         };
      };
}

#endif // NONDIMESIONAL_HPP
