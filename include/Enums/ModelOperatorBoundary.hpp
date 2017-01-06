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

namespace QuICC {

      /**
       * @brief Simple struct holding the mapping for the model operators boundary types
       */
      struct ModelOperatorBoundary {
         /**
          * @brief Enums of different model operator boundary types
          */
         enum Id {
            /// Solver operator with a boundary condition (Tau or Galerkin)
            SOLVER_HAS_BC = 0,
            /// Solver operator without Tau boundary but with Galerkin condition (ie RHS part)
            SOLVER_NO_TAU,
            /// Special setup for Galerkin stencil matrix
            STENCIL,
            /// Setup for operators projecting from field values to solver RHS (ie, QI, explicit linear terms)
            FIELD_TO_RHS,
         };
      };
}

#endif // MODELOPERATORBOUNDARY_HPP
