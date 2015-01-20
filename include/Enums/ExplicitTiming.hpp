/**
 * @file ExplicitTiming.hpp
 * @brief Definition of some useful enums for the explicit linear timing with respect to nonlinear calculation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EXPLICITTIMING_HPP
#define EXPLICITTIMING_HPP

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
    * @brief Some timing enums for equation explicit linear terms
    */
   struct ExplicitTiming {

      /**
      * @name Enum for the timining of the equation with respect to prognostic solve
      */
      enum Id {
         /// Use explicit linear matrix on field values
         LINEAR = 0,
         /// Use explicit linear matrix on nonlinear values (ie. as nondiagonal QI terms)
         NONLINEAR,
      };
   };
}

#endif // EXPLICITTIMING_HPP
