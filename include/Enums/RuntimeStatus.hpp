/**
 * @file RuntimeStatus.hpp
 * @brief Definition of some useful enums for runtime status 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef RUNTIMESTATUS_HPP
#define RUNTIMESTATUS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @brief Simple struck holding the runtime statuses
    */
   struct RuntimeStatus {

      /**
       * @name Enum for useful runtime statuses
       */
      enum Id {
         /// Keep going
         GOON = 0, 
         /// Abort simulation
         STOP = 2
      };
   };
}

#endif // RUNTIMESTATUS_HPP
