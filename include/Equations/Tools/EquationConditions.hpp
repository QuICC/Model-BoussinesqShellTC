/** 
 * @file EquationConditions.hpp
 * @brief Implementation of equation condition functors
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONCONDITIONS_HPP
#define EQUATIONCONDITIONS_HPP

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

namespace Conditions {

   /**
    * @brief Condition functor for prognostic equation
    */
   struct IsPrognostic
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for diagnostic equation
    */
   struct IsDiagnostic
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for trivial equation
    */
   struct IsTrivial
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

   /**
    * @brief Condition functor for wrapper
    */
   struct IsWrapper
   {
      /**
       * @brief Check scalar equation
       *
       * @param eqA Scalar equation
       */
      bool operator()(SharedIScalarEquation eqA);

      /**
       * @brief Check vector equation
       *
       * @param eqA vector equation
       */
      bool operator()(SharedIVectorEquation eqA);
   };

}
}
}

#endif // EQUATIONCONDITIONS_HPP
