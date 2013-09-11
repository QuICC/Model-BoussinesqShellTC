/**
 * @file RuntimeStatus.hpp
 * @brief Definition of some useful enums for runtime status 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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

namespace GeoMHDiSCC {

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
         STOP
      };
   };
}

#endif // RUNTIMESTATUS_HPP
