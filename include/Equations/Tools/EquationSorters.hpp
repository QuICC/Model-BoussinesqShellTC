/** 
 * @file EquationTools.hpp
 * @brief Implementation of equation sorting functors
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONSORTERS_HPP
#define EQUATIONSORTERS_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Sorters {

   /**
    * @brief Sorting functor by equation type
    */
   class EquationType
   {
      public:
         /**
          * @brief Sort scalar equations by equation type
          *
          * @param eqA Left scalar equation
          * @param eqB Right scalar equation
          */
         bool operator()(SharedIScalarEquation eqA, SharedIScalarEquation eqB);

         /**
          * @brief Sort vector equations by equation type 
          *
          * @param eqA Left vector equation
          * @param eqB Right vector equation
          */
         bool operator()(SharedIVectorEquation eqA, SharedIVectorEquation eqB);

      private:
         /**
          * @brief Convert enum equation type to integer for scalar equation
          *
          * @param eqA Scalar equation
          */
         int computeEquationType(SharedIScalarEquation eqA);

         /**
          * @brief Convert enum equation type to integer for vector equation
          *
          * @param eqA Vector equation
          */
         int computeEquationType(SharedIVectorEquation eqA);

   };

}
}
}

#endif // EQUATIONSORTERS_HPP
