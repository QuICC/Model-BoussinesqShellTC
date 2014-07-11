/**
 * @file SolveTiming.hpp
 * @brief Definition of some useful enums for the solve timing with respect to prognostic solve 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SOLVETIMING_HPP
#define SOLVETIMING_HPP

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
    * @brief Some timing enums for equation solves
    */
   struct SolveTiming {

      /**
      * @name Enum for the timining of the equation with respect to prognostic solve
      */
      enum Id {
         /// Solve before timestep
         BEFORE = 0,
         /// Solve is a prognostic equation
         PROGNOSTIC,
         /// Solve after
         AFTER,
      };
   };
}

#endif // SOLVETIMING_HPP
