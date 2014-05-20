/**
 * @file ModelOperator.hpp
 * @brief Definition of some useful enums for the model operators 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MODELOPERATOR_HPP
#define MODELOPERATOR_HPP

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
       * @brief Simple struct holding the mapping for the model operators
       */
      struct ModelOperator {
         /**
          * @brief Enums of the different model operators
          */
         enum Id {
            /// Quasi inverse operator
            QI = 0,
            /// Time operator (block diagonal)
            TIME,
            /// Implicitly timestepped linear operator
            IMPLICIT_LINEAR,
            /// Explicitly timestepped linear operator
            EXPLICIT_LINEAR,
         };
      };
}

#endif // MODELOPERATOR_HPP
