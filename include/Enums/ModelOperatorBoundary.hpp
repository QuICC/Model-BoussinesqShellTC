/**
 * @file ModelOperatorBoundary.hpp
 * @brief Definition of some useful enums for the model operators boundary types
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MODELOPERATORBOUNDARY_HPP
#define MODELOPERATORBOUNDARY_HPP

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
       * @brief Simple struct holding the mapping for the model operators boundary types
       */
      struct ModelOperatorBoundary {
         /**
          * @brief Enums of different model operator boundary types
          */
         enum Id {
            /// Operator has a boundary condition (Tau or Galerkin)
            HAS_BC = 0,
            /// Operator has no Tau boundary but Galerkin should be applied
            NO_TAU,
            /// Operator is boundary free (no Tau and no Galerkin)
            NO_BC,
         };
      };
}

#endif // MODELOPERATORBOUNDARY_HPP
