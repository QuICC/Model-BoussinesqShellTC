/**
 * @file CartesianExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANEXACTSTATEIDS_HPP
#define CARTESIANEXACTSTATEIDS_HPP

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
   struct CartesianExactStateIds
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
         // TTT states
         POLYPOLYPOLY = 10,  // Polynomial, Polynomial, Polynomial
         // TFT states
         POLYCOSPOLY = 20, // Polynomial, Cosine, Polynomial
         POLYSINPOLY,      // Polynomial, Sine, Polynomial
         // TFF states
         POLYCOSCOS = 30,  // Polynomial, Cosine, Cosine
         POLYSINSIN,       // Polynomial, Sine, Sine
         POLYSINCOS,       // Polynomial, Sine, Cosine
         POLYCOSSIN,       // Polynomial, Cosine, Sine
         // FFF states
         COSCOSCOS = 50,   // Cosine, Cosine, Cosine
         SINSINSIN,        // Sine, Sine, Sine
         COSCOSSIN,        // Cosine, Cosine, Sine
         SINSINCOS,        // Sine, Sine, Cosine
         COSSINSIN,        // Cosine, Sine, Sine
         SINCOSCOS,        // Sine, Cosine, Cosine
         COSSINCOS,        // Cosine, Sine, Cosine
         SINCOSSIN,        // Sine, Cosine, Sine
         // Special
         SPECIAL1 = 1000,  // Kind of a place holder for tests
         SPECIAL2,         // Kind of a place holder for tests
         SPECIAL3,         // Kind of a place holder for tests
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
   };

}
}

#endif // CARTESIANEXACTSTATEIDS_HPP
