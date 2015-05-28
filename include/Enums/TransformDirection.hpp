/**
 * @file TransformDirection.hpp
 * @brief Definition of some useful enums for the direction of the transforms needed for the nonlinear calculations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMDIRECTION_HPP
#define TRANSFORMDIRECTION_HPP

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
    * @brief Transform direction list
    */
   struct TransformDirection {

      /**
      * @name Enum for the direction of the transform
      */
      enum Id {
         /// Forward direction of transform: physical -> spectral
         FORWARD = 0,
         /// Backward direction of transform: spectral -> physical
         BACKWARD,
      };
   };
}

#endif // TRANSFORMDIRECTION_HPP
