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
         enum Id {EKMAN, ROBERTS, RAYLEIGH, ROSSBY, MAGEKMAN, PRANDTL, MAGPRANDTL, CHI, GAMMA};
      };
}

#endif // NONDIMESIONAL_HPP
