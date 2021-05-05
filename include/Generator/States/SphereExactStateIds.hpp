/**
 * @file SphereExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHEREEXACTSTATEIDS_HPP
#define SPHEREEXACTSTATEIDS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

namespace Equations {

   /**
    * @brief Class to hold the list of possible exact states
    */
   struct SphereExactStateIds
   {
      /// Polynomial approximation to Cosine
      static const MHDFloat PCOS;

      /// Polynomial approximation to Sine
      static const MHDFloat PSIN;

      /**
       * @brief Enums for the avaialable exact states
       */
      enum Id {
         NOTUSED = -1,  // Initialisation state (do NOT use)
         // Special states
         CONSTANT = 0,  // All constant

         HARMONIC = 10, // Generate spherical harmonic state

         BENCHTEMPC1 = 50, // Initial temperature perturbation state for full sphere benchmark C1
         BENCHVELC1,       // Initial velocity perturbation state for full sphere benchmark C1
         BENCHVELC2,       // Initial velocity perturbation state for full sphere benchmark C2
	 VALIDATION_ENSTROPHY,	// State designed to validate the enstrophy writer
         BENCHMAGC2,       // Initial magnetic perturbation state for full sphere benchmark C2
      };

      /**
       * @brief Compute spherical harmonic physical values
       */
      static Array sph_harmonic(const MHDComplex amplitude, const int l, const int m, const MHDFloat theta, const Array& phi);
   };

}
}

#endif // SPHEREEXACTSTATEIDS_HPP
