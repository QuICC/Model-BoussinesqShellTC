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

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Class to hold the list of possible exact states
    */
   struct SphereExactStateIds
   {
      /// Polynomial approximation to Cosine
      static const MHDFloat PCOS = 99999;

      /// Polynomial approximation to Sine
      static const MHDFloat PSIN = -99999;

      /**
       * @brief Enums for the avaialable exact states
       */
      enum Id {
         // Special states
         CONSTANT = 0,  // All constant
      };
   };

}
}

#endif // SPHEREEXACTSTATEIDS_HPP