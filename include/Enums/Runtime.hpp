/** \file Runtime.hpp
 *  \brief Definition of some useful enums for runtime information
 *
 *  \mhdBug Needs test
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
       *
       *  - GOON: Keep going
       *  - STOP: Abort simulation
       */
      enum Id {GOON = 0, STOP};
   };
}
}

#endif // RUNTIME_HPP
