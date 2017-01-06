/**
 * @file CylinderExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDEREXACTSTATEIDS_HPP
#define CYLINDEREXACTSTATEIDS_HPP

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
   struct CylinderExactStateIds
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
         POLYCOSPOLY = 10, // Polynomial, Cosine, Polynomial
         POLYSINPOLY,      // Polynomial, Sine, Polynomial
         SPEC_UNIT, // Spectral space solution: unit spectrum of max 5 modes
      };

      /**
       * @brief Compute even periodic mode
       */
      static MHDFloat cos(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta);

      /**
       * @brief Compute odd periodic mode
       */
      static MHDFloat sin(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat theta);

      /**
       * @brief Compute polynomial mode
       */
      static MHDFloat poly(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x);

      /**
       * @brief Compute Chebyshev mode
       */
      static MHDFloat chebyshev(const MHDFloat amplitude, const MHDFloat mode, const MHDFloat x);
   };

}
}

#endif // CYLINDEREXACTSTATEIDS_HPP
