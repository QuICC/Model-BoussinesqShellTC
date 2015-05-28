/**
 * @file ShellExactStateIds.hpp
 * @brief Class to hold the list of possible exact states
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHELLEXACTSTATEIDS_HPP
#define SHELLEXACTSTATEIDS_HPP

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
   struct ShellExactStateIds
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

         HARMONIC = 10, // Generate spherical harmonic state

         TOROIDAL = 20, // Toroidal Y_l^m
         POLOIDAL,      // Poloidal Y_l^m

         TORPOL = 30, // Toroidal+Poloidal state

         BENCHTEMPC1 = 50, // Initial temperature perturbation state for Christensen's benchmark C1
         BENCHVELC1,       // Initial velocity perturbation state for Christensen's benchmark C1
         BENCHMAGC1,       // Initial magnetic perturbation state for Christensen's benchmark C1
      };

      /**
       * @brief Compute spherical harmonic physical values
       */
      static Array sph_harmonic(const MHDComplex amplitude, const int l, const int m, const MHDFloat theta, const Array& phi);
   };

}
}

#endif // SHELLEXACTSTATEIDS_HPP
