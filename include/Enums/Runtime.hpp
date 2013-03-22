/** \file Runtime.hpp
 *  \brief Definition of some useful enums for runtime information
 */

#ifndef RUNTIME_HPP
#define RUNTIME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Runtime {

   /**
    * @brief Simple struck holding the runtime statuses
    */
   struct Status {

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
}

#endif // RUNTIME_HPP
