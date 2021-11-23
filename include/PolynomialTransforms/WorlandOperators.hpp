/**
 * @file WorlandOperators.hpp
 * @brief Implementation of some Worland operators
 */

#ifndef WORLANDOPERATORS_HPP
#define WORLANDOPERATORS_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

namespace Polynomial {

   /**
    * @brief Implementation of some Worland operators
    */
   class WorlandOperators
   {
      public:
         /**
          * @brief Compute \f$W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void integralR3(Matrix& op, const int n, const int l);

      private:
         /**
          * @brief Constructor
          */
         WorlandOperators() = default;

         /**
          * @brief Destructor
          */
         ~WorlandOperators() = default;

   };
}
}

#endif // WORLANDOPERATORS_HPP
